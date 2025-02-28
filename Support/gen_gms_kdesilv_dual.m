% Generate Gaussian Mixture Samples from a Dual GMM using Kernel Density 
% Estimation and Silverman's Rule of Thumb
function [x_out, P_out, w_out] = gen_gms_kdesilv_dual(model, w1, m1, P1, w2, m2, P2, num_par)
%
    W1 = sum(w1);
    W2 = sum(w2);
    Z = W1 + W2;

    w1cumsum = cumsum(w1/sum(w1));
    w2cumsum = cumsum(w2/sum(w2));

    if num_par > 0
        x_out = zeros(model.x_dim,num_par);
        for ii = 1:num_par
            u = rand(model.stream, [1,1]);
            if u < W1/Z % select from first mixture
                v = rand(model.stream, [1,1]);
                idx = find(w1cumsum >= v, 1, 'first');
                mk_ii = m1(:,idx);
                sqrt_Pk_ii = sqrtm(P1(:,:,idx));
            else % select from second mixture
                v = rand(model.stream, [1,1]);
                idx = find(w2cumsum >= v, 1, 'first');
                mk_ii = m2(:,idx);
                sqrt_Pk_ii = sqrtm(P2(:,:,idx));
            end
            x_out(:,ii)   = mk_ii + sqrt_Pk_ii * randn(model.stream, [model.x_dim,1]);
        end 
        w_out = Z * ones(num_par,1) ./ num_par; % uniform weights
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