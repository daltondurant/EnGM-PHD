% Generate Gaussian Mixture Samples using Kernel Density Estimation and
% Silverman's Rule of Thumb
function [x_out, P_out, w_out] = gen_gms_kdesilv(model,w,m,P,num_par)
%
    if num_par > 0
        u = rand(model.stream,[num_par,1]);
        wcumsum = cumsum(w/sum(w));
        x_out = zeros(model.x_dim,num_par);
        P_out = zeros(model.x_dim,model.x_dim,num_par);
        for ii = 1:num_par
            idx = find(wcumsum >= u(ii), 1, 'first');
            sqrt_Pk = sqrtm(P(:,:,idx));
            x_out(:,ii)   = m(:,idx) + sqrt_Pk * randn(model.stream, [model.x_dim,1]) ;
            P_out(:,:,ii) = P(:,:,idx); % if sampling from set, then return component's covariance
        end
        w_out = sum(w) * ones(num_par,1) ./ num_par; % uniform weights
        betaS_scale = 1; % Silverman's rule of thumb heuristic scaling coefficient (>1 more conservative / <1 more confident)
        if num_par > 1
            betaS = (betaS_scale/ceil(sum(w_out)))*(4/(num_par*(model.x_dim+2))) ^ (2/(model.x_dim+4)); % Silverman's rule of thumb
            mux = mean(x_out,2); ex = x_out - mux; P = betaS * ((ex * ex') / (num_par-1) ); % Pbar (sample covariance)
            P_out = repmat(P,1,1,num_par);
        else
            P_out = [];
        end
    else
        x_out = [];
        P_out = [];
        w_out = [];
    end
%
end