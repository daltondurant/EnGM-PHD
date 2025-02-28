% Gaussian merging
function [w_new,x_new,P_new]= gaus_merge(w,x,P,varargin)
%
    merge_threshold = 4; % DEFAULT 4-sigma merging threshold

    %--- loading any optional arguments
    while ~isempty(varargin)
        switch lower(varargin{1}) % switch to lowercase only
            case 'mergethreshold'
                  merge_threshold = varargin{2};
            otherwise
                  error(['Unexpected option: ' varargin{1}])
        end
        varargin(1:2) = [];
    end

    L= length(w); 
    x_dim= size(x,1);
    I= 1:L;
    el= 1;
    
    if all(w==0)
        w_new = [];
        x_new = [];
        P_new = [];
        return;
    end
    
    while ~isempty(I)
        [~,j]= max(w); 
        j = j(1);
        Ij= []; 
        iPt= inv(P(:,:,j));
        w_new(el,1)= 0; 
        x_new(:,el)= zeros(x_dim,1); 
        P_new(:,:,el)= zeros(x_dim,x_dim);
        for i= I
            val= (x(:,i)-x(:,j))'*iPt*(x(:,i)-x(:,j));
            if val <= merge_threshold
                Ij= [ Ij i ];
            end
        end
    
        w_new(el,1)   = sum(w(Ij));
        wmat = repmat(w(Ij)',[x_dim,1]);
        x_new(:,el)  = sum(wmat.*x(:,Ij),2);
        wtemp = reshape(w(Ij),[1,1,size(w(Ij))]);
        wmat = repmat(wtemp,[x_dim,x_dim,1]);
        P_new(:,:,el) = sum(wmat.*P(:,:,Ij),3);
        
        x_new(:,el)= x_new(:,el)/w_new(el);
        P_new(:,:,el)= P_new(:,:,el)/w_new(el);


        I= setdiff(I,Ij);
        w(Ij)= -1;
        el= el+1;
    end
%
end
