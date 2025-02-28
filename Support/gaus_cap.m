% Gaussian capping
function [w_new,x_new,P_new]= gaus_cap(w,x,P,varargin)
%
    L_max = 100; % DEFAULT limit on number of Gaussians

    %--- loading any optional arguments
    while ~isempty(varargin)
        switch lower(varargin{1}) % switch to lowercase only
            case 'lmax'
                  L_max = varargin{2};
            otherwise
                  error(['Unexpected option: ' varargin{1}])
        end
        varargin(1:2) = [];
    end

    if length(w) > L_max
        [~,idx]= sort(w,1,'descend');
        w_new= w(idx(1:L_max)); 
        w_new = w_new * (sum(w)/sum(w_new));
        x_new= x(:,idx(1:L_max));
        P_new= P(:,:,idx(1:L_max));
    else
        x_new = x;
        P_new = P;
        w_new = w;
    end
%
end

