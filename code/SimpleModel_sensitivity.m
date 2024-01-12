tic
warning('off', 'MATLAB:polyshape:repairedBySimplify');
warning('off', 'MATLAB:polyshape:boundary3Points');
%% Parameters
alpha_1_min = 0.05; % equilibrium ratio of internal P to dissolved+adsorbed P  [-] %0.125
alpha_1_max = 0.25;
alpha_2_min = 0.002; % fraction of dissolved+adsorbed pool that is dissolved at equilibrium [-] % 0.005
alpha_2_max = 0.02;
alpha_3_min = 0.1; % rate of internal diffusion[1/T] % 0.25
alpha_3_max = 0.3;
Pup_min = 20; % used to set alpha_4 parameter, average uptake [kg/ha] for initial condition % 27 if subsoil supply via fert
Pup_max = 25; % 32 if subsoil supply via fert
% alpha_4_min = 33; % maximum plant uptake of dissolved P [P] % 36
% alpha_4_max = 36; % above 35.5 recurrent sets sometimes not in safe seet ~36.25 gives high rate of such
alpha_5_min = 20; % for critical P level = 58 kg/ha; Michaelis-Menten parameter for plant uptake of dissolved P [P] % 40
alpha_5_max = 45; % for critical P level = 135 kg/ha
% alpha_5_min = -10; % Michaelis-Menten parameter for plant uptake of dissolved P [P] % 40
% alpha_5_max = 10;
DRP_PP_min = 0.3; % Ratio of dissolved p to particulate p edge-of-field losses
DRP_PP_max = 1.2;
% beta_1_min = 0.025;  % fraction of dissolved P lost with runoff per unit flow [-] % 0.15 or 0.175
% beta_1_max = 0.25;
% beta_2_min = 0.000025; % fraction of particulate P lost with runoff per unit flow [-] (for 1 m depth, Q_bar=4, 0.00025 -> 1.0 mm loss) % 0.00015 or 0.000175
% beta_2_max = 0.00025;
Loss_min = 1.6; % total edge-of-field losses [kg/ha]
Loss_max = 2.2;
% Q_bar_min = -20; % mean flow % 4
% Q_bar_max = 20;
num_sets = 10000;
params = lhsdesign(num_sets,7);
% params = lhsdesign(num_sets,8); % OLD
u_min = 0; % "natural" P weathering rate % 7 if incorporating subsoil supply here
u_max = 30; % max 30 kg/ha fertilizer input % +7 if incorporating subsoil supply here
x1_init = 200; % 226 for 20 cm, 200 for 17.75 cm (7 inches)
% ave_load = 2.0; % 1.923 if using Q_bar=4 beta_1 = 0175 and 20 cm, 2.0 if using 17.75 cm and 0.175
% ave_yield = 30;

%% Run code and identify hitting times for all parameterizations
hit_vector_opt = zeros(num_sets,1);
hit_vector_stquo = zeros(num_sets,1);
hit_vector_constrained = zeros(num_sets,1);
delta_time = zeros(num_sets,1);
Yr45Reach = repmat(polyshape, 1, num_sets);
horizon = 45; % number of seasons or years to simulate
f = waitbar(0, 'Starting scenarios');
%%%%%debug begin 
% for zz=44:44
%%%%% debug end
for zz=1:num_sets
    % Parameter value, tested for sensitivity
    alpha_1 = alpha_1_min + (params(zz,1))*(alpha_1_max-alpha_1_min);
    x2_init = x1_init/alpha_1 + 42; % a little above equilibrium since fertilizer applications are lower that decades previous
    alpha_2 = alpha_2_min + (params(zz,2))*(alpha_2_max-alpha_2_min);
    alpha_3 = alpha_3_min + (params(zz,3))*(alpha_3_max-alpha_3_min);
    alpha_5 = alpha_5_min + (params(zz,5))*(alpha_5_max-alpha_5_min); 
    Pup = Pup_min + (params(zz,4))*(Pup_max-Pup_min); 
    alpha_4 = Pup *(x1_init+alpha_5)/x1_init;
%     alpha_4 = alpha_4_min + (params(zz,4))*(alpha_4_max-alpha_4_min);
%     alpha_5 = (alpha_4 - ave_yield)*(x1_init/ave_yield) * (1+0.01*(params(zz,5)-0.5)*(alpha_5_max-alpha_5_min));  % 30 is approximate P uptake in data
    DRP_PP = DRP_PP_min + (params(zz,6))*(DRP_PP_max-DRP_PP_min);
%     B1B2 = ((1+ (1/alpha_1)) / alpha_2) * DRP_PP;
    B1B2 = ((1-alpha_2)*x1_init+x2_init)/(alpha_2*x1_init) * DRP_PP;
    beta_1 = 0.175;
    beta_2 = beta_1/B1B2;
%     beta_1 = beta_1_min + (params(zz,6))*(beta_1_max-beta_1_min);
%     beta_2 = beta_2_min + (params(zz,7))*(beta_2_max-beta_2_min);
    Loss = Loss_min + (params(zz,7))*(Loss_max-Loss_min);
    Q_bar = Loss/((beta_1*alpha_2+beta_2*(1-alpha_2))*x1_init + beta_2*x2_init);  
%     Q_bar = ave_load/((beta_1*alpha_2+beta_2*(1-alpha_2))*x1_init + beta_2*x2_init) * (1+0.01*(params(zz,8)-0.5)*(Q_bar_max-Q_bar_min));  % 1.8 is approximate annual load in data
    % Q_bar = Q_bar_min + (params(zz,8))*(Q_bar_max-Q_bar_min);

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
        run('SimpleModel_switchregions_sens.m'); % identify u_max region
        run('SimpleModel_hittime.m'); % calculate hitting time for no intermediate contstraint
        if hit_time==0
            stop=0;
            mm = 8;
            while stop==0
                Singular_array_all = subtract(Singular_array_all,Singular_array(mm));
                mm = mm-1;
                run('SimpleModel_hittime.m'); % calculate hitting time for no intermediate contstraint
                if (hit_time>0) || (mm==3)
                    stop=1;
                end
            end
        end
        hit_vector_opt(zz) = hit_time;
        run('SimpleModel_stquotime.m'); % calculate hitting time for no intermediate contstraint
        hit_vector_stquo(zz) = hit_time;
        run('SimpleModel_constrainedtime.m'); % calculate hitting time for no intermediate contstraint
        hit_vector_constrained(zz) = hit_time;
    else
        hit_vector_opt(zz) = NaN;
        hit_vector_stquo(zz) = NaN;
        hit_vector_constrained(zz) = NaN;
    end
    waitbar(zz/num_sets, f, sprintf('Progress: %d %%', floor(zz/num_sets*100)));
end
close(f)
delta_time = hit_vector_constrained - hit_vector_opt;
%%
toc

%% Plots
figure; hist(hit_vector_opt(find(hit_vector_opt>0)),(10:3:130))
title('', 'fontsize', 26)
xlabel('Hitting time (years)', 'fontsize', 26)
ylabel('Count', 'fontsize', 26)
ax = gca; % axes handle
ax.FontSize = 18; 
figure; hist(hit_vector_constrained(find(hit_vector_constrained>0)),50)
title('', 'fontsize', 26)
xlabel('Hitting time (years)', 'fontsize', 26)
ylabel('Count', 'fontsize', 26)
ax = gca; % axes handle
ax.FontSize = 18; 
figure; hist(delta_time(find(hit_vector_constrained>0)),250)
title('', 'fontsize', 26)
xlabel({'Difference in optimal versus','constrained hitting time (years)'}, 'fontsize', 26)
ylabel('Count', 'fontsize', 26)
ax = gca; % axes handle
ax.FontSize = 18; 
xlim([0 100]); % to make relevant portion more visible
figure; hist(hit_vector_stquo(find(hit_vector_stquo>0)),50)
title('', 'fontsize', 26)
xlabel('Hitting time (years)', 'fontsize', 26)
ylabel('Count', 'fontsize', 26)
ax = gca; % axes handle
ax.FontSize = 18; 
figure; scatter(params(:,1),hit_vector_opt)
title('\alpha_1', 'fontsize', 26)
xlabel('Normalized parameter value', 'fontsize', 26)
ylabel('Minimum hitting time', 'fontsize', 26)
ax = gca; % axes handle
ax.FontSize = 18; 
figure; scatter(params(:,2),hit_vector_opt)
title('\alpha_2', 'fontsize', 26)
xlabel('Normalized parameter value', 'fontsize', 26)
ylabel('Minimum hitting time', 'fontsize', 26)
ax = gca; % axes handle
ax.FontSize = 18; 
figure; scatter(params(:,3),hit_vector_opt)
title('\alpha_3', 'fontsize', 26)
xlabel('Normalized parameter value', 'fontsize', 26)
ylabel('Minimum hitting time', 'fontsize', 26)
ax = gca; % axes handle
ax.FontSize = 18; 
figure; scatter(params(:,4),hit_vector_opt)
title('Average annual plant P uptake', 'fontsize', 26)
xlabel('Normalized parameter value', 'fontsize', 26)
ylabel('Minimum hitting time', 'fontsize', 26)
ax = gca; % axes handle
ax.FontSize = 18; 
figure; scatter(params(:,5),hit_vector_opt)
title('\alpha_5', 'fontsize', 26)
xlabel('Normalized parameter value', 'fontsize', 26)
ylabel('Minimum hitting time', 'fontsize', 26)
ax = gca; % axes handle
ax.FontSize = 18; 
figure; scatter(params(:,6),hit_vector_opt)
title('Ratio of dissolve:particulate P losses', 'fontsize', 26)
xlabel('Normalized parameter value', 'fontsize', 26)
ylabel('Minimum hitting time', 'fontsize', 26)
ax = gca; % axes handle
ax.FontSize = 18; 
figure; scatter(params(:,7),hit_vector_opt)
title('Average annual edge-of-field P losses', 'fontsize', 26)
xlabel('Normalized parameter value', 'fontsize', 26)
ylabel('Minimum hitting time', 'fontsize', 26)
ax = gca; % axes handle
ax.FontSize = 18; 
figure; scatter(params(:,1),delta_time)
title('\alpha_1', 'fontsize', 26)
xlabel('Normalized parameter value', 'fontsize', 26)
ylabel({'Difference in optimal versus','constrained hitting time (years)'}, 'fontsize', 26)
ax = gca; % axes handle
ax.FontSize = 18; 
ylim([0 max(delta_time)]); % to make relevant portion more visible
figure; scatter(params(:,5),delta_time)
title('\alpha_5', 'fontsize', 26)
xlabel('Normalized parameter value', 'fontsize', 26)
ylabel({'Difference in optimal versus','constrained hitting time (years)'}, 'fontsize', 26)
ax = gca; % axes handle
ax.FontSize = 18; 
ylim([0 max(delta_time)]); % to make relevant portion more visible

%% Statistics
nan_opt = length(find(isnan(hit_vector_opt)));
zero_opt = length(find(hit_vector_opt==0));
mean_opt = mean(hit_vector_opt(~isnan(hit_vector_opt)));
med_opt = median(hit_vector_opt(~isnan(hit_vector_opt)));
min_opt = min(hit_vector_opt(~isnan(hit_vector_opt)));
max_opt = max(hit_vector_opt(~isnan(hit_vector_opt)));

nan_constr = length(find(isnan(hit_vector_constrained)));
zero_constr = length(find(hit_vector_constrained==0));
mean_constr = mean(hit_vector_constrained(~isnan(hit_vector_constrained)));
med_constr = median(hit_vector_constrained(~isnan(hit_vector_constrained)));
min_constr = min(hit_vector_constrained(~isnan(hit_vector_constrained)));
max_constr = max(hit_vector_constrained(~isnan(hit_vector_constrained)));

nan_stquo = length(find(isnan(hit_vector_stquo)));
zero_stquo = length(find(hit_vector_stquo==0));
mean_stquo = mean(hit_vector_stquo(find(hit_vector_stquo>0)));
med_stquo = median(hit_vector_stquo(find(hit_vector_stquo>0)));
min_stquo = min(hit_vector_stquo(find(hit_vector_stquo>0)));
max_stquo = max(hit_vector_stquo(find(hit_vector_stquo>0)));

mean_delta = mean(delta_time(~isnan(hit_vector_opt)));
med_delta = median(delta_time(~isnan(hit_vector_opt)));
min_delta = min(delta_time(~isnan(hit_vector_opt)));
max_delta = max(delta_time(~isnan(hit_vector_opt)));

