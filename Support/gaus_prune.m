% Gaussian pruning
function [w_new,x_new,P_new]= gaus_prune(w,x,P,varargin)
%
    elim_threshold = 1e-5; % DEFAULT pruning threshold

    %--- loading any optional arguments
    while ~isempty(varargin)
        switch lower(varargin{1}) % switch to lowercase only
            case 'elimthreshold'
                  elim_threshold = varargin{2};
            otherwise
                  error(['Unexpected option: ' varargin{1}])
        end
        varargin(1:2) = [];
    end

    idx= find( w > elim_threshold );
    w_new= w(idx);
    x_new= x(:,idx);
    P_new= P(:,:,idx);
%
end

