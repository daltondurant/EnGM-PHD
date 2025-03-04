%% engmphd
%**************************************************************************
% author: dalton durant
% email: ddurant@utexas.edu, thedaltondurant@gmail.com
%
% USE AT YOU OWN DISCRETION
%
% This code is free to use, all I ask is you cite me when reasonable.
%**************************************************************************
% The Ensemble Gaussian Mixture Probability Hypothesis Density filter.
%
% SOURCE: 
% [1] Dalton Durant and Renato Zanetti. "Kernel-Based Ensemble Gaussian 
% Mixture Probability Hypothesis Density Filter." In 28th International 
% Conference on Information Fusion, Rio de Janeiro, Brazil. 2025.
%
function [model,meas,est] = run_filter(stream,cfig,model,meas,est)
%
    % --- filter parameters
    run_flag= 'disp';   % 'disp' or 'silence' for on the fly output
    J_rsp   = 250;      % resampling fixed number of particles
    J_birth = 10;       % number of births per time step 
        
    % --- output variables
    est.X = cell(meas.K,1);  % states
    est.N = zeros(meas.K,1); % cardinality 
    est.G = zeros(meas.K,1); % number of Gaussians
    
    model.name  = 'EnGM-PHD';
    model.stream = stream;

    % --- initial prior
    w_predict = 1e-16;
    m_predict = model.m0 + sqrtm(model.P0)*randn(stream,[model.x_dim,1]);
    P_predict = model.P0;

    betaS_scale = 1; % Silverman's rule of thumb heuristic scaling coefficient (>1 more conservative / <1 more confident)

    % --- recursive filtering
    for k = 1:meas.K % time loop
        % 1. Time handling
        model = cfig.time(model, k);

        % 2. Births
        %
        % directly bakes birth model into prior intensity GMM
        [m_predict, P_predict, w_predict] = gen_gms_kdesilv_dual(model, w_predict, m_predict, P_predict, ...
                                                model.w_birth, model.m_birth, model.P_birth, J_rsp + J_birth);
        %
        %{
        % constructs a GMM from birth model first, then bakes into prior intensity GMM
        [m_birth, P_birth, w_birth] = gen_gms_kdesilv(model,model.w_birth,model.m_birth,model.P_birth,J_birth); 
        [m_predict, P_predict, w_predict] = gen_gms_kdesilv_dual(model, w_predict, m_predict, P_predict, ...
                                                w_birth, m_birth, P_birth, J_rsp + J_birth);
        %}
                                                                                                              
        % 3. Gating
        if model.gate_flag
            Zk = gate_meas(cfig, model, meas.Z{k}, m_predict, P_predict);        
        else
            Zk = meas.Z{k};
        end
            
        % 4. Update
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

        % 5. Resampling
        [m_update, ~, w_update] = gen_gms_kdesilv(model,w_update,m_update,P_update,J_rsp);

        % 6. Predict
        % --- survivors
        m_predict = cfig.fprop(model, m_update, 'noiseless');
        betaS = (betaS_scale/ceil(sum(w_update)))*(4/(J_rsp*(model.x_dim+2))) ^ (2/(model.x_dim+4)); % Silverman's rule of thumb
        mum = mean(m_predict,2); ex = m_predict - mum; Pum = ((ex * ex') / (J_rsp-1)); % Pbar (sample covariance)
        P_predict = repmat((betaS*Pum)+model.Q,1,1,J_rsp);
        w_predict= model.P_S*w_update; % surviving weights

        % 7. State Extraction   
        if sum(w_update) > 0.5
            [x_c,I_c]= our_kmeans(model.stream,m_update,w_update,1);
            est.N(k)= 0;
            for j=1:size(x_c,2)
                if sum(w_update(I_c{j})) > .5
                    est.X{k}= [ est.X{k} x_c(:,j) ];
                    est.N(k)= est.N(k)+1;
                end
            end
        else
            est.N(k)= 0; est.X{k}= [];
        end
        est.G(k) = length(w_update);

        % 8. Diagnostics
        if ~strcmp(run_flag,'silence')
            disp([' time= ',num2str(k),...
             ' #est mean=' num2str(sum(w_update),4),... % estimated mean number of targets
             ' #est card=' num2str(est.N(k),4) ]);      % estimated cardinality
        end
    end % end time loop
%
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Helper functions
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