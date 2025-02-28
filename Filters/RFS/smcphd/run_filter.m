%% smcphd
%**************************************************************************
% author: dalton durant
% email: ddurant@utexas.edu, thedaltondurant@gmail.com
%
% USE AT YOU OWN DISCRETION
%
% This code is free to use, all I ask is you cite me when reasonable.
%**************************************************************************
% The Sequential Monte Carlo Probability Hypothesis Density filter with random
% resampling.
%
% SOURCE: 
% [1] Ba-Ngu Vo, Sumeetpal Singh, and Arnaud Doucet. "Sequential Monte 
% Carlo methods for Bayesian multi-target filtering with random finite 
% sets." IEEE Trans. Aerospace and Electronic Systems. 2005.
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
    
    model.name  = 'SMC-PHD';
    model.stream = stream;

    % --- initial prior
    w_update = 1e-16;
    m_update = model.m0 + sqrtm(model.P0)*randn(stream,[model.x_dim,1]);

    % --- recursive filtering
    for k = 1:meas.K % time loop
        % 1. Time handling
        model = cfig.time(model, k);

        % 2. Predict
        % --- survivors
        m_predict = predict(cfig,model,m_update); % surviving components
        w_predict= model.P_S*w_update; % surviving weights
        
        % 3. Births
        [m_birth, ~, w_birth] = gen_gms(model,model.w_birth,model.m_birth,model.P_birth,J_birth); 
        % --- append
        m_predict = cat(2,m_predict,m_birth); 
        w_predict = cat(1,w_predict,w_birth); 
            
        % 4. Update
        Zk = meas.Z{k};
        m = size(Zk,2); % number of measurements
        % --- missed detections
        w_update = model.Q_D*w_predict;
        m_update = m_predict;
        % --- BPF update
        if m~=0
            for ii=1:m
                U_temp = compute_likelihood(cfig,model,Zk(:,ii),model.R,m_predict);
                
                % this is a numerically stable way to do this using the log-sum-exp trick
                log_w_temp = log(model.P_D) + log(w_predict(:)) + U_temp;
                log_sum_w_temp = max(log_w_temp) + log(sum(exp(log_w_temp - max(log_w_temp))));
                log_denom_w_temp = log(model.lambda_c*model.pdf_c + exp(log_sum_w_temp));
                log_w_temp = log_w_temp - log_denom_w_temp;
                w_update = w_update + exp(log_w_temp);
            end
        end
        
        % 5. Resampling
        idx = randsample(model.stream,length(w_update),J_rsp,true,w_update); 
        w_update = sum(w_update)*ones(J_rsp,1)/J_rsp;
        m_update = m_update(:,idx); % this is a naive way of doing this

        % 6. State Extraction   
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
        
        % 7. Diagnostics
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
% Predict States with Random Noise from Process Noise Model
function Xkp1 = predict(cfig,model,Xk)   
% 
    Xkp1 = zeros(size(Xk));
    for idxp=1:size(Xk,2)
        Xkp1(:,idxp) = cfig.fprop(model, Xk(:,idxp), 'noise');
    end
%
end

% Compute Likelihood
function U_update = compute_likelihood(cfig,model,y,R,m)
%      
    det_R= prod(eig(R)); iR = pinv(R); % precompute
    U_update = zeros(size(m,2),1);
    for jj=1:size(m,2)
        ybar = cfig.h(model,m(:,jj),'noiseless');
        U_update(jj) = -(y-ybar)' * (iR * (y-ybar)) / 2 ...
                       - log(det_R) / 2 ...
                       - log(2*pi) * size(y,1) / 2;
    end
%
end