% Gaussian splitting
function [w_new,x_new,P_new]= gaus_split(w,x,P,varargin)
%
    split_displacement = 0.5; % DEFAULT split displacement

    %--- loading any optional arguments
    while ~isempty(varargin)
        switch lower(varargin{1}) % switch to lowercase only
            case 'splitdisplacement'
                  split_displacement = varargin{2};
            otherwise
                  error(['Unexpected option: ' varargin{1}])
        end
        varargin(1:2) = [];
    end
    
    if all(w==0)
        x_new = []; P_new = []; w_new = [];
        return;
    end
    
    x_new = []; P_new = []; w_new = [];
    for ii = 1:length(w)
        nu = split_displacement; % displacement
        mu = x(:,ii);
        Sigma = P(:,:,ii);
        [V,D,~] = eig(Sigma); 
        [lambda, max_idx] = max(diag(D));
        u = V(:,max_idx);
    
        % 3 component split
        mu1 = mu + nu*sqrt(lambda)*u;
        mu2 = mu - nu*sqrt(lambda)*u;
        mu3 = mu;
        Sigma1 = Sigma - (nu^2)*lambda*(u*u')/3;
        Sigma2 = Sigma - (nu^2)*lambda*(u*u')/3;
        Sigma3 = Sigma - (nu^2)*lambda*(u*u')/3;
        w1 = 1/6 * w(ii);
        w2 = 1/6 * w(ii);
        w3 = 4/6 * w(ii);

        x_new = cat(2, x_new, mu1);
        x_new = cat(2, x_new, mu2);
        x_new = cat(2, x_new, mu3);
        P_new = cat(3, P_new, Sigma1);
        P_new = cat(3, P_new, Sigma2);
        P_new = cat(3, P_new, Sigma3);
        w_new = cat(1, w_new, w1);
        w_new = cat(1, w_new, w2);
        w_new = cat(1, w_new, w3);
    end
end
