%% PARAMETERS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Quick changing useful parameters
model.dt           = 1;  % sampling rate
truth.K = 100;                        % total time steps 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Basic Parameters

% --- state space and measurement space
model.x_dim= 6;   %dimension of state vector
model.z_dim= 3;   %dimension of observation vector

% --- dynamics
model.T0 = model.dt; % sim/filter rate
time_array = (0:1:truth.K) * model.T0;
model.len_time = length(time_array);
model.Time = time_array'; 
proc_u = 1; % process noise uncertainty 
model.B= proc_u*eye(model.x_dim);         
model.Q= model.B*model.B';  % process noise covariance

% --- measurements
model.rIs = zeros(3,1); % pos of sensor
model.vIs = zeros(3,1); % vel of sensor
rho_u = 1; % range uncertainty 
angle_u = 0.5*pi/180; % angle uncertainity [rad]
az_u = angle_u; % azimuth uncertainty [rad] 
el_u = angle_u; % elevation uncertainty [rad] 
model.D= diag([rho_u; az_u; el_u]);
model.R= model.D*model.D'; % covariance for observation noise 


% --- clutter
model.lambda_c = 10; % poisson average rate of uniform clutter per scan
model.range_c= [ 0 200; 0 200; 0 400];           
model.pdf_c= 1/prod(model.range_c(:,2)-model.range_c(:,1)); %uniform clutter density


% --- births
model.birth_pos_u = 50;       
model.birth_vel_u = 5; 
model.w_birth= 1/100; 
model.B_birth= diag([model.birth_pos_u; model.birth_pos_u; model.birth_pos_u; ...
                     model.birth_vel_u; model.birth_vel_u; model.birth_vel_u]);
%model.P_birth = model.Q; % everyone wins scenario
model.P_birth= model.B_birth*model.B_birth'; % hard case
model.m_birth= [ 75; 75; 150; 0; 0; 0];


% --- gating
model.gate_flag = 0;   % gating on or off 1/0
model.P_G = 0.99;      % gate size in percentage
model.gamma_G = chi2inv(model.P_G,model.z_dim);  % gate inv chi^2 dn gamma value 


% --- detection parameters
model.P_D= .98;   % probability of detection in measurements
model.Q_D= 1-model.P_D; % probability of missed detection in measurements


% --- survival parameters
model.P_S= .99; % probability of survival
model.Q_S= 1-model.P_S; % probability of death


% --- initial conditions
model.m0 = [0;0;0;0;0;0];
model.P0 = eye(model.x_dim);


% --- truth and initial states
truth.total_tracks = 2;    % number of targets
% target initial conditions
X0(:,1) = [50, 50, 50, 1/2, 1/2, 2];
X0(:,2) = [100, 100, 50, -1/2, -1/2, 2];
X = X0;
truth.X = cell(model.len_time,1);
model.T = model.T0;
for tk = 1:model.len_time        
    % store data
    truth.X{tk} = X;
    % propagate true chaser and target states
    if tk < model.len_time
        X = cfig.fprop(model, X, 'noiseless');
    end
end




