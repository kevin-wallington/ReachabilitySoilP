warning('off', 'MATLAB:polyshape:repairedBySimplify');
warning('off', 'MATLAB:polyshape:boundary3Points');
%% Parameters
horizon = 100; % number of seasons or years to simulate
alpha_1 = 0.125; % equilibrium ratio of internal P to dissolved+adsorbed P  [-] 
alpha_2 = 0.005; % fraction of dissolved+adsorbed pool that is dissolved at equilibrium [-] 
alpha_3 = 0.25; % rate of internal diffusion[1/T] 
alpha_4 = 27; % maximum plant uptake of dissolved P [P] 
alpha_5 = 35; % Michaelis-Menten parameter for plant uptake of dissolved P [P]
beta_1 = 0.175; % fraction of dissolved P lost with runoff per unit flow [-] % 0.15
beta_2 = 0.000175; % fraction of particulate P lost with runoff per unit flow [-] (for 1 m depth, Q_bar=4, 0.00025 -> 1.0 mm of sediment loss) % 0.00015
Q_bar = 4; % mean flow 
x1_init = 200; % [kg/ha] % 226
x2_init = x1_init/alpha_1 + 42; % a little above equilibrium since fertilizer applications are lower that decades previous

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
u_min = 0; % "natural" P weathering rate... actually subsoil supply rate??? 7 for subsoil supply, 1 for weathering % 7 if including subsoil supply
u_max = 30; % max 30 kg/ha fertilizer input - +7 for subsoil supply, +1 for weathering % +7 if including subsoil supply

%% Lower limit for x1 in state space
mask_x = 30;

%% Compute Original Safe Set
Def_lim = 0.75; % nutrient deficit limit (i.e. 75% of max)
Load_lim = 0.55*Q_bar*((beta_1*alpha_2+beta_2*(1-alpha_2))*x1_init + beta_2*x2_init); % 226 and 1850 are estimated IC
x1_min = Def_lim/(1-Def_lim)*alpha_5;
x2_max = (Load_lim - Q_bar*(beta_1*alpha_2+beta_2*(1-alpha_2))*x1_min) / (Q_bar*beta_2);
x1_max = Load_lim /(Q_bar*(beta_1*alpha_2+beta_2*(1-alpha_2)));
cx = [x1_min, x1_min, x1_max]; % corner points, x-coord
cy = [0, x2_max, 0]; % corner points, y-coord
poly_points = [cx; cy]';
pgon=polyshape(poly_points);
SafeSet_orig=pgon;

%% Compute Safe Set and Reachable Set with Interventions
beta_1 = 0.175; % fraction of dissolved P lost with runoff per unit flow [-] % 0.15
beta_2 = 0.000125; % fraction of particulate P lost with runoff per unit flow [-] (for 1 m depth, Q_bar=4, 0.00025 -> 1.0 mm of sediment loss) % 0.00015
x1_min = Def_lim/(1-Def_lim)*alpha_5;
x2_max = (Load_lim - Q_bar*(beta_1*alpha_2+beta_2*(1-alpha_2))*x1_min) / (Q_bar*beta_2);
x1_max = Load_lim /(Q_bar*(beta_1*alpha_2+beta_2*(1-alpha_2)));
cx = [x1_min, x1_min, x1_max]; % corner points, x-coord
cy = [0, x2_max, 0]; % corner points, y-coord
poly_points = [cx; cy]';
pgon=polyshape(poly_points);
SafeSet_b2=pgon;
SafeSet = pgon; 
run('SimpleModel_targetonly.m'); % don't need recurent set, safe set from above
run('SimpleModel_reachability.m'); % complete
Target_b2 = SafeInv;
Reach40_b2 = Reach_array(40);
Reach18_b2 = Reach_array(18);
Reach10_b2 = Reach_array(10);
Reach3_b2 = Reach_array(3);
Reach2_b2 = Reach_array(2);
Reach1_b2 = Reach_array(1);

beta_1 = 0.125; % fraction of dissolved P lost with runoff per unit flow [-] % 0.15
beta_2 = 0.000175; % fraction of particulate P lost with runoff per unit flow [-] (for 1 m depth, Q_bar=4, 0.00025 -> 1.0 mm of sediment loss) % 0.00015
x1_min = Def_lim/(1-Def_lim)*alpha_5;
x2_max = (Load_lim - Q_bar*(beta_1*alpha_2+beta_2*(1-alpha_2))*x1_min) / (Q_bar*beta_2);
x1_max = Load_lim /(Q_bar*(beta_1*alpha_2+beta_2*(1-alpha_2)));
cx = [x1_min, x1_min, x1_max]; % corner points, x-coord
cy = [0, x2_max, 0]; % corner points, y-coord
poly_points = [cx; cy]';
pgon=polyshape(poly_points);
SafeSet_b1=pgon;
SafeSet = pgon; 
run('SimpleModel_targetonly.m'); % don't need recurent set, safe set from above
run('SimpleModel_reachability.m'); % complete
Target_b1 = SafeInv;
Reach40_b1 = Reach_array(40);
Reach33_b1 = Reach_array(33);
Reach10_b1 = Reach_array(10);
Reach3_b1 = Reach_array(3);
Reach2_b1 = Reach_array(2);
Reach1_b1 = Reach_array(1);

figure
hold on
plot(SafeSet_b2,'FaceColor',[0.85 0.85 0.85], 'FaceAlpha',1);
plot(SafeSet_b1,'FaceColor',[0.15 0.15 0.15], 'FaceAlpha',0.7);
plot(SafeSet_orig,'FaceColor',[0 0 0], 'FaceAlpha',1);
ax = gca; % axes handle
ax.FontSize = 18; 
title('', 'fontsize', 26)
xlabel('Labile P (kg/ha)', 'fontsize', 26)
ylabel('Stable P (kg/ha)', 'fontsize', 26)
lgd = legend('\beta_2 reduction','\beta_1 reduction','No reduction', 'Location', 'Northeast', 'fontsize', 26);
lgd.Color = 'white';
legend boxoff
xlim([100 400])
ylim([0 1500])

figure
hold on
plot(Reach40_b2,'FaceColor',[0.8 0.8 1.0],'EdgeColor',[1 1 1],'FaceAlpha',1);
plot(Reach18_b2,'FaceColor',[0.7 0.7 0.9],'EdgeColor',[1 1 1],'FaceAlpha',1);
plot(Reach10_b2,'FaceColor',[0.4 0.4 0.7],'EdgeColor',[1 1 1],'FaceAlpha',1);
plot(Reach3_b2,'FaceColor',[0.2 0.2 0.5],'EdgeColor',[1 1 1],'FaceAlpha',1);
plot(Reach2_b2,'FaceColor',[0.1 0.1 0.3],'EdgeColor',[1 1 1],'FaceAlpha',1);
plot(Reach1_b2,'FaceColor',[0 0.0 0.1],'EdgeColor',[1 1 1],'FaceAlpha',1);
plot(Target_b2,'FaceColor',[0.5 0.5 0.5],'EdgeColor',[1 1 1],'FaceAlpha',1);
% plot initial condition estimate
plot([x1_init],[x2_init],'pentagram','MarkerFaceColor',[0.75 0 0.15],'MarkerEdgeColor','white', 'MarkerSize',50,'linewidth',2);
% plot formatting
xlim([mask_x 220])
ylim([mask_y 1800])
ax = gca; % axes handle
ax.FontSize = 18; 
title('\beta_2', 'fontsize', 26)
xlabel('Labile P (kg/ha)', 'fontsize', 26)
ylabel('Stable P (kg/ha)', 'fontsize', 26)
lgd = legend('Reachable in 40 years','18 years','10 years','3 years','2 years','1 year', 'Target', 'Location', 'Northwest', 'fontsize', 26);
lgd.Color = 'white';
legend boxon

figure
hold on
plot(Reach40_b1,'FaceColor',[0.8 0.8 1.0],'EdgeColor',[1 1 1],'FaceAlpha',1);
plot(Reach33_b1,'FaceColor',[0.7 0.7 0.9],'EdgeColor',[1 1 1],'FaceAlpha',1);
plot(Reach10_b1,'FaceColor',[0.4 0.4 0.7],'EdgeColor',[1 1 1],'FaceAlpha',1);
plot(Reach3_b1,'FaceColor',[0.2 0.2 0.5],'EdgeColor',[1 1 1],'FaceAlpha',1);
plot(Reach2_b1,'FaceColor',[0.1 0.1 0.3],'EdgeColor',[1 1 1],'FaceAlpha',1);
plot(Reach1_b1,'FaceColor',[0 0.0 0.1],'EdgeColor',[1 1 1],'FaceAlpha',1);
plot(Target_b1,'FaceColor',[0.5 0.5 0.5],'EdgeColor',[1 1 1],'FaceAlpha',1);
% plot initial condition estimate
plot([x1_init],[x2_init],'pentagram','MarkerFaceColor',[0.75 0 0.15],'MarkerEdgeColor','white', 'MarkerSize',50,'linewidth',2);
% plot formatting
xlim([mask_x 220])
ylim([mask_y 1800])
ax = gca; % axes handle
ax.FontSize = 18; 
title('\beta_1', 'fontsize', 26)
xlabel('Labile P (kg/ha)', 'fontsize', 26)
ylabel('Stable P (kg/ha)', 'fontsize', 26)
lgd = legend('Reachable in 40 years','33 years','10 years','3 years','2 years','1 year', 'Target', 'Location', 'Northwest', 'fontsize', 26);
lgd.Color = 'white';
legend boxon