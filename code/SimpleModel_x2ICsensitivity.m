tic
warning('off', 'MATLAB:polyshape:repairedBySimplify');
warning('off', 'MATLAB:polyshape:boundary3Points');
%% Parameters
horizon = 100; % number of seasons or years to simulate
alpha_1 = 0.125; % equilibrium ratio of internal P to dissolved+adsorbed P  [-] 
alpha_2 = 0.005; % fraction of dissolved+adsorbed pool that is dissolved at equilibrium [-] 
alpha_3 = 0.25; % rate of internal diffusion[1/T] 
alpha_4 = 36; % maximum plant uptake of dissolved P [P] 
alpha_5 = 40; % Michaelis-Menten parameter for plant uptake of dissolved P [P]
beta_1 = 0.175; % fraction of dissolved P lost with runoff per unit flow [-] % 0.15
beta_2 = 0.000175; % fraction of particulate P lost with runoff per unit flow [-] (for 1 m depth, Q_bar=4, 0.00025 -> 1.0 mm of sediment loss) % 0.00015
Q_bar = 4; % mean flow 
DRP_PP = 0.6;
Loss = 2.0; % [kg/ha] average annual losses currently
x1_init = 200; % 
x2_init_vec = (1000:50:5000); % 81 ICs, at 2550 (yields alpha_1=0.0784), x2 change only allows min fert to hit target with no switch
num_sets = length(x2_init_vec);
u_min = 7; % "natural" P weathering rate
u_max = 30 + 7; % max 30 kg/ha fertilizer input

%% Run code and identify hitting times for all parameterizations
hit_vector_opt_ic = zeros(num_sets,1);
hit_vector_stquo_ic = zeros(num_sets,1);
hit_vector_constrained_ic = zeros(num_sets,1);
hit_vector_opt_alpha = zeros(num_sets,1);
hit_vector_stquo_alpha = zeros(num_sets,1);
hit_vector_constrained_alpha = zeros(num_sets,1);
% x1_plot = zeros(num_sets,1);
% x2_plot = zeros(num_sets,1);
Yr45Reach = repmat(polyshape, 1, num_sets);
horizon = 45; % number of seasons or years to simulate
f = waitbar(0, 'Starting scenarios');
%%%%%debug begin 
% for zz=14:14
%%%%% debug end
for zz=1:num_sets
    x2_init = x2_init_vec(zz);
    Q_bar = Loss/((beta_1*alpha_2+beta_2*(1-alpha_2))*x1_init + beta_2*x2_init);
    
    %%% For just changing the initial x2 (and beta_2)
    alpha_1 = 0.125;
    B1B2 = ((1-alpha_2)*x1_init + x2_init)/(alpha_2*x1_init) * DRP_PP;
    beta_1 = 0.175;
    beta_2 = beta_1/B1B2;
    Q_bar = Loss/((beta_1*alpha_2+beta_2*(1-alpha_2))*x1_init + beta_2*x2_init);
    % Dynamics (linear)
    A = zeros(2,2); 
    A_loc = zeros(2,2); 
    A(1,1) = 1 - alpha_3 - (beta_1*alpha_2 + beta_2*(1-alpha_2))*Q_bar;
    A(1,2) = alpha_3*alpha_1;
    A(2,1) = alpha_3;
    A(2,2) = 1 - alpha_3*alpha_1 - beta_2*Q_bar;
    
    % For output generation
    C = zeros(2,2); % state -> output, deterministic
    C(1,1) = 1; % adsorbed+dissolved state
    C(2,1) = (beta_1*alpha_2 + beta_2*(1 - alpha_2))*Q_bar; % losses from x1
    C(2,2) = beta_2*Q_bar; % losses from x2
    
    run('SimpleModel_checkSafeRecur.m'); 
    if ~isempty(x1i)
        run('SimpleModel_targetonly.m'); % don't need recurent set, safe set from above
        run('SimpleModel_invlimit.m'); % to set mask x
        run('SimpleModel_reachability_sens.m'); % necessary to identify singular regions - omits plots, and doesn't use umaxback
        Singular_array_all = subtract(Singular_array_all,Singular_array(8)); % these get wonky for standard params
        Singular_array_all = subtract(Singular_array_all,Singular_array(7));
        run('SimpleModel_switchregions_sens.m'); % identify u_max region
        run('SimpleModel_hittime.m'); % calculate hitting time for no intermediate contstraint
        hit_vector_opt_ic(zz) = hit_time;
        run('SimpleModel_stquotime.m'); % calculate hitting time for no intermediate contstraint
        hit_vector_stquo_ic(zz) = hit_time;
        run('SimpleModel_constrainedtime.m'); % calculate hitting time for no intermediate contstraint
        hit_vector_constrained_ic(zz) = hit_time;
    else
        hit_vector_opt_ic(zz) = NaN;
        hit_vector_stquo_ic(zz) = NaN;
        hit_vector_constrained_ic(zz) = NaN;
    end 
    
    %%% For changing the initial x2 (and beta_2) AND alpha_1
    alpha_1 = x1_init/x2_init;
    B1B2 = ((1+ (1/alpha_1)) / alpha_2) * DRP_PP;
    beta_1 = 0.175;
    beta_2 = beta_1/B1B2;
    Q_bar = Loss/((beta_1*alpha_2+beta_2*(1-alpha_2))*x1_init + beta_2*x2_init);
    % Dynamics (linear)
    A = zeros(2,2); 
    A_loc = zeros(2,2); 
    A(1,1) = 1 - alpha_3 - (beta_1*alpha_2 + beta_2*(1-alpha_2))*Q_bar;
    A(1,2) = alpha_3*alpha_1;
    A(2,1) = alpha_3;
    A(2,2) = 1 - alpha_3*alpha_1 - beta_2*Q_bar;
    
    % For output generation
    C = zeros(2,2); % state -> output, deterministic
    C(1,1) = 1; % adsorbed+dissolved state
    C(2,1) = (beta_1*alpha_2 + beta_2*(1 - alpha_2))*Q_bar; % losses from x1
    C(2,2) = beta_2*Q_bar; % losses from x2
    
    run('SimpleModel_checkSafeRecur.m'); 
    if ~isempty(x1i)
        run('SimpleModel_targetonly.m'); % don't need recurent set, safe set from above
        run('SimpleModel_invlimit.m'); % to set mask x
        run('SimpleModel_reachability_sens.m'); % necessary to identify singular regions - omits plots, and doesn't use umaxback
        Singular_array_all = subtract(Singular_array_all,Singular_array(8)); % these get wonky for standard params
        Singular_array_all = subtract(Singular_array_all,Singular_array(7));
        run('SimpleModel_switchregions_sens.m'); % identify u_max region
        run('SimpleModel_hittime.m'); % calculate hitting time for no intermediate contstraint
        hit_vector_opt_alpha(zz) = hit_time;
        run('SimpleModel_stquotime.m'); % calculate hitting time for no intermediate contstraint
        hit_vector_stquo_alpha(zz) = hit_time;
        run('SimpleModel_constrainedtime.m'); % calculate hitting time for no intermediate contstraint
        hit_vector_constrained_alpha(zz) = hit_time;
    else
        hit_vector_opt_alpha(zz) = NaN;
        hit_vector_stquo_alpha(zz) = NaN;
        hit_vector_constrained_alpha(zz) = NaN;
    end 
    
    waitbar(zz/num_sets, f, sprintf('Progress: %d %%', floor(zz/num_sets*100)));
end
close(f)
%%
toc

%% Plots
figure
hold on
scatter(x2_init_vec, hit_vector_opt_ic, 'filled', 'MarkerFaceColor','r');
scatter(x2_init_vec, hit_vector_opt_alpha, 'filled', 'MarkerFaceColor','b');
