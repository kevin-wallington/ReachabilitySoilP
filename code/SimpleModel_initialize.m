% Initialize parameters and matrices for dynamics
% Note, this uses 2D model, see retired versions for simulating 3D model
% (3rd state is in-stream P)
%% Parameters
% Note: parameters for in-stream processing not verified
horizon = 100; % number of seasons or years to simulate
alpha_1 = 0.125; % equilibrium ratio of internal P to dissolved+adsorbed P  [-] 
alpha_2 = 0.005; % fraction of dissolved+adsorbed pool that is dissolved at equilibrium [-] 
alpha_3 = 0.25; % rate of internal diffusion[1/T] 
alpha_4 = 27; % 35; % 36 % maximum plant uptake of dissolved P [P] % 35 and 35 gives average annual yield ~ 29.7 kg/ha ; 27 for new framing of subsoil supply
alpha_5 = 35; % 40 % Michaelis-Menten parameter for plant uptake of dissolved P [P], 35 gives critical bray P1 = 27.2 mg/kg
beta_1 = 0.175; % fraction of dissolved P lost with runoff per unit flow [-] % 0.15
beta_2 = 0.000175; % fraction of particulate P lost with runoff per unit flow [-] (for 1 m depth, Q_bar=4, 0.00025 -> 1.0 mm of sediment loss) % 0.00015
Q_bar = 4; % mean flow 
x1_init = 200; % [kg/ha] % 226
x2_init = x1_init/alpha_1 + 42; % a little above equilibrium since fertilizer applications are lower that decades previous
% ave_load = 1.8; %[kg/ha] annually
% Q_bar = ave_load/((beta_1*alpha_2+beta_2*(1-alpha_2))*x1_init + beta_2*x2_init);
%% Dynamics (linear, deterministic only)
A = zeros(2,2); 
A_loc = zeros(2,2); 
A(1,1) = 1 - alpha_3 - (beta_1*alpha_2 + beta_2*(1-alpha_2))*Q_bar;
A(1,2) = alpha_3*alpha_1;
A(2,1) = alpha_3;
A(2,2) = 1 - alpha_3*alpha_1 - beta_2*Q_bar;

%% For output generation
C = zeros(2,2); % state -> output, deterministic
C(1,1) = 1; % adsorbed+dissolved state
C(2,1) = (beta_1*alpha_2 + beta_2*(1 - alpha_2))*Q_bar; % losses from x1
C(2,2) = beta_2*Q_bar; % losses from x2

%% Control vector
u_min = 0;%7; % "natural" P weathering rate... actually subsoil supply rate??? 7 for subsoil supply, 1 for weathering; 0 for new framing of subsoil supply
u_max = 30;% + 7; % max 30 kg/ha fertilizer input - +7 for subsoil supply, +1 for weathering

%% Lower limit for x1 in state space
mask_x = 30;

