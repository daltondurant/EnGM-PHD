% Generate Gaussian Mixture Samples
function [x_out, P_out, w_out] = gen_gms(model,w,m,P,num_par)
%
    if num_par > 0
        u = rand(model.stream,[num_par,1]);
        wcumsum = cumsum(w/sum(w));
        x_out = zeros(model.x_dim,num_par);
        P_out = zeros(model.x_dim,model.x_dim,num_par);
        for ii = 1:num_par
            idx = find(wcumsum >= u(ii), 1, 'first');
            x_out(:,ii)   = m(:,idx) + sqrtm(P(:,:,idx)) * randn(model.stream, [model.x_dim,1]) ;
            P_out(:,:,ii) = P(:,:,idx); % this is a naive way of doing this
        end
        w_out = sum(w) * ones(num_par,1) ./ num_par;
    else
        x_out = [];
        P_out = [];
        w_out = [];
    end

%
end
