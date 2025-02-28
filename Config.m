classdef Config
    properties
        % constants
        options
        nMonte

    end
    methods
        function obj = Config()
            % creates obj object
        end

        function model = time(obj, model, k)
            if k == 1
                model.T = 0;
            else
                model.T = model.T0;
            end
        end
        
        % measurement function
        function y = h(obj, model, X, eta)
            if ~isnumeric(eta)
                if strcmp(eta,'noise')
                    eta = model.D*randn(model.stream,[size(model.D,2),size(X,2)]);
                elseif strcmp(eta,'noiseless')
                    eta = zeros(size(model.D,2),size(X,2));
                end
            end
            
            if isempty(X)
                y= [];
            else
                y = zeros(model.z_dim, size(X,2));

                for i = 1:size(X,2)
                    rI = X(1:3,i);
                    rho = rI-model.rIs;
                    rhonorm = sqrt(sum(rho.^2));                   % range
                    az   = atan2(rho(2,:),rho(1,:));               % azimuth [rad]
                    el   = asin(rho(3,:)/rhonorm);                 % elevation [rad]
                    y(:,i) = [rhonorm; az; el] + eta(:,i);
                end
            end
        end
        
        % measurement jacobian
        function H = H(obj, model, X)
            %{
            syms rIs1 rIs2 rIs3 vIs1 vIs2 vIs3 x1 x2 x3 xdot1 xdot2 xdot3 real
            rI = [x1; x2; x3];
            vI = [xdot1; xdot2; xdot3];
            rIs = [rIs1; rIs2; rIs3];
            vIs = [vIs1; vIs2; vIs3];
            rho = rI - rIs;
            rhonorm = sqrt(sum(rho.^2)); % range
            rhonormdot = (rho)' * (vI - vIs) / rhonorm; % range rate
            alpha   = atan2(rho(2,:),rho(1,:)); % right ascension [rad]
            delta   = asin(rho(3)/rhonorm); % declination [rad]
            %y = [rhonorm; alpha; delta];
            %y = [rhonorm; rhonormdot; alpha; delta];
            y = [alpha; delta];
            jacobian(y,[x1,x2,x3,xdot1,xdot2,xdot3])
            %}
            x1 = X(1,:);
            x2 = X(2,:);
            x3 = X(3,:);
            xdot1 = X(4,:);
            xdot2 = X(5,:);
            xdot3 = X(6,:);
            r = model.rIs;
            v = model.vIs;
            rIs1 = r(1);
            rIs2 = r(2);
            rIs3 = r(3);
            vIs1 = v(1);
            vIs2 = v(2);
            vIs3 = v(3);
            % range, angles
            H = [-(0.5000*(2*rIs1 - 2*x1))/((rIs1 - x1)^2 + (rIs2 - x2)^2 + (rIs3 - x3)^2)^0.5000,                                                                                          -(0.5000*(2*rIs2 - 2*x2))/((rIs1 - x1)^2 + (rIs2 - x2)^2 + (rIs3 - x3)^2)^0.5000,                                                                                                                                                    -(0.5000*(2*rIs3 - 2*x3))/((rIs1 - x1)^2 + (rIs2 - x2)^2 + (rIs3 - x3)^2)^0.5000, 0, 0, 0;
                  (rIs2 - x2)/((rIs1 - x1)^2 + (rIs2 - x2)^2),                                                                                                                              -(rIs1 - x1)/((rIs1 - x1)^2 + (rIs2 - x2)^2),                                                                                                                                                                                                                                   0, 0, 0, 0;
                 -(0.5000*(2*rIs1 - 2*x1)*(rIs3 - x3))/((1 - (rIs3 - x3)^2/((rIs1 - x1)^2 + (rIs2 - x2)^2 + (rIs3 - x3)^2))^0.5000*((rIs1 - x1)^2 + (rIs2 - x2)^2 + (rIs3 - x3)^2)^1.5000), -(0.5000*(2*rIs2 - 2*x2)*(rIs3 - x3))/((1 - (rIs3 - x3)^2/((rIs1 - x1)^2 + (rIs2 - x2)^2 + (rIs3 - x3)^2))^0.5000*((rIs1 - x1)^2 + (rIs2 - x2)^2 + (rIs3 - x3)^2)^1.5000), (1/((rIs1 - x1)^2 + (rIs2 - x2)^2 + (rIs3 - x3)^2)^0.5000 - (0.5000*(2*rIs3 - 2*x3)*(rIs3 - x3))/((rIs1 - x1)^2 + (rIs2 - x2)^2 + (rIs3 - x3)^2)^1.5000)/(1 - (rIs3 - x3)^2/((rIs1 - x1)^2 + (rIs2 - x2)^2 + (rIs3 - x3)^2))^0.5000, 0, 0, 0];
        end

        % dynamics function
        function f = f(obj, ~, X, model)
            x   =X(1,:);
            y   =X(2,:);
            z   =X(3,:);
            xdot=X(4,:);
            ydot=X(5,:);
            zdot=X(6,:);
            
            f = zeros(size(X,1),size(X,2));
            f(1,:)=xdot;
            f(2,:)=ydot;
            f(3,:)=zdot;
            f(4,:)=0;
            f(5,:)=0;
            f(6,:)=0;

        end

        % dynamics jacobian 
        function F = F(obj, X, model)
            F = [zeros(3,3) eye(3);
                zeros(3,6)];
        end
        
        % time propagation
        function Xkp1 = fprop(obj, model, Xk, V)  
            
            if ~isnumeric(V)
                if strcmp(V,'noise')
                    V= model.B*randn(size(model.B,2),size(Xk,2));
                elseif strcmp(V,'noiseless')
                    V = zeros(size(model.B,2),size(Xk,2));
                end
            end
            
            if isempty(Xk)
                Xkp1 = [];
            else 
                [~,m] = eDP54(@(T,X) obj.f(T, X, model), [0 model.T], Xk, obj.options);
                Xkp1 = m + V;
                
            end

        end
        
       
    end % end methods
end % end class