%% gmphd
%**************************************************************************
% author: dalton durant
% email: ddurant@utexas.edu, thedaltondurant@gmail.com
%
% USE AT YOU OWN DISCRETION
%
% This code is free to use, all I ask is you cite me when reasonable.
%**************************************************************************
% The Gaussian Mixture Probability Hypothesis Density filter with pruning,
% merging, and capping.
%
% SOURCE: 
% [1] Ba-Ngu Vo and Wing-Kin Ma. "The Gaussian mixture probability
% hypothesis density filter." IEEE Transactions on Signal Processing. 2006.
% 
function [model,meas,est] = run_filter(stream,cfig,model,meas,est)
%
    % --- filter parameters
    run_flag= 'disp';             % 'disp' or 'silence' for on the fly output
    J_birth = 10;                 % number of births per time step 
    elim_threshold  = 1e-5;       % pruning threshold
    merge_threshold = 4;          % merging threshold
    cap_limit = 250;              % capping limit
        
    % --- output variables
    est.X = cell(meas.K,1);  % states
    est.N = zeros(meas.K,1); % cardinality 
    est.G = zeros(meas.K,1); % number of Gaussians
    
    model.name  = 'GM-PHD';
    model.stream = stream;

    % --- initial prior
    w_update = 1e-16;
    m_update = model.m0 + sqrtm(model.P0)*randn(stream,[model.x_dim,1]);
    P_update = model.P0;

    % --- recursive filtering
    for k = 1:meas.K % time loop
        % 1. Time handling
        model = cfig.time(model, k);

        % 2. Predict
        % --- survivors
        [m_predict,P_predict] = predict(cfig,model,m_update,P_update); % surviving components
        w_predict= model.P_S*w_update; % surviving weights
        
        % 3. Births
        [m_birth, P_birth, w_birth] = gen_gms(model,model.w_birth,model.m_birth,model.P_birth,J_birth); 
        % --- append
        m_predict = cat(2,m_predict,m_birth); 
        P_predict = cat(3,P_predict,P_birth);
        w_predict = cat(1,w_predict,w_birth);  
                                                                                                     
        % 4. Gating
        if model.gate_flag
            Zk = gate_meas(cfig, model, meas.Z{k}, m_predict, P_predict);        
        else
            Zk = meas.Z{k};
        end
            
        % 5. Update
        m = size(Zk,2); % number of measurements
        % --- missed detections
        w_update = model.Q_D*w_predict;
        m_update = m_predict;
        P_update = P_predict;    
        % --- KF update and compute hypotheses
        if m~=0
            for ii=1:m
                [m_temp,P_temp,U_temp] = update(cfig,model,Zk(:,ii),model.R,m_predict,P_predict);
                
                % this is a numerically stable way to do this using the log-sum-exp trick
                log_w_temp = log(model.P_D) + log(w_predict(:)) + U_temp;
                log_sum_w_temp = max(log_w_temp) + log(sum(exp(log_w_temp - max(log_w_temp))));
                log_denom_w_temp = log(model.lambda_c*model.pdf_c + exp(log_sum_w_temp));
                log_w_temp = log_w_temp - log_denom_w_temp;
                w_temp = exp(log_w_temp);

                w_update = cat(1,w_update,w_temp);
                m_update = cat(2,m_update,m_temp);
                P_update = cat(3,P_update,P_temp);
            end
        end

        % 6. Pruning, merging, and capping
        if m~=0
            [w_update,m_update,P_update]= gaus_prune(w_update,m_update,P_update,...
                'ElimThreshold',elim_threshold);
            [w_update,m_update,P_update]= gaus_merge(w_update,m_update,P_update,...
                'MergeThreshold',merge_threshold);   
            [w_update,m_update,P_update]= gaus_cap(w_update,m_update,P_update,...
                'Lmax',cap_limit);
        end

        % 7. State extraction [1] (not used by filter)
        idx= find(w_update > 0.5 );
        for ii = 1:length(idx)
            repeat_num_targets= round(w_update(idx(ii)));
            est.X{k} = [ est.X{k} repmat(m_update(:,idx(ii)),[1,repeat_num_targets]) ];
            est.N(k) = est.N(k) + repeat_num_targets;
        end
        est.G(k) = length(w_update);
        
        % 8. Diagnostics
        if ~strcmp(run_flag,'silence')
            disp([' time= ',num2str(k),...
             ' #est mean=' num2str(sum(w_update),4),... % estimated mean number of targets
             ' #est card=' num2str(est.N(k),4)]);       % estimated cardinality 
        end
    end % end time loop
%
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Helper functions
% Predict States and Covariances for the Extended Kalman Filter (EKF)
function [Xkp1,Pkp1] = predict(cfig,model,Xk,Pk)   
% 
    Xkp1 = zeros(size(Xk));
    Pkp1 = zeros(size(Pk));
    for idxp=1:size(Xk,2)
        Xkp1(:,idxp) = cfig.fprop(model, Xk(:,idxp), 'noiseless');
        Phi = expm(cfig.F(Xk(:,idxp), model)*model.T); % approximation for STM for small time steps
        %-- add scaled noise 
        Pkp1(:,:,idxp) = Phi*Pk(:,:,idxp)*Phi' + model.Q;
        Pkp1(:,:,idxp) = (Pkp1(:,:,idxp) + Pkp1(:,:,idxp)')/2;
    end
%
end

% Update States and Covariances for the Extended Kalman Filter (EKF)
function [m_update,P_update,U_update] = update(cfig,model,y,R,m,P)
%      
    U_update = zeros(size(m,2),1);
    m_update = zeros(size(m,1),size(m,2));
    P_update = zeros(size(m,1),size(m,1),size(m,2));
    for jj=1:size(m,2)
        % individual EKF update
        ybar = cfig.h(model,m(:,jj),'noiseless');
        Hj   = cfig.H(model,m(:,jj));
        Pxxj = P(:,:,jj);
        Pxyj = Pxxj * Hj';
        Pyyj= Hj*Pxxj*Hj' + R; 
        Pyyj= (Pyyj + Pyyj')/2;   % additional step to avoid numerical problems
        det_Pyyj = prod(eig(Pyyj)); iPyyj = pinv(Pyyj);
        Kj = Pxyj * iPyyj;
        m_update(:,jj)  = m(:,jj) + Kj*(y-ybar);
        Ij = eye(size(m,1));
        P_update(:,:,jj) = (Ij - Kj*Hj) * Pxxj * (Ij - Kj*Hj)' + Kj * R * Kj'; % Joseph form
        
        % weight update
        U_update(jj) = -(y-ybar)' * (iPyyj * (y-ybar)) / 2 ...
                       - log(det_Pyyj) / 2 ...
                       - log(2*pi) * size(y,1) / 2;
    end
%
end