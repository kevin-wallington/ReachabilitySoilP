tic
warning('off', 'MATLAB:polyshape:repairedBySimplify');
warning('off', 'MATLAB:polyshape:boundary3Points');
%% Set up target based on status quo
run('SimpleModel_initialize.m'); % complete
Def_lim = 0.75; % nutrient deficit limit (i.e. 75% of max)
Load_lim = 0.55*Q_bar*((beta_1*alpha_2+beta_2*(1-alpha_2))*x1_init + beta_2*x2_init); % 226 and 1850 are estimated IC

%% Parameters
num_sets = 20;
beta_1_min = 0.125; % 0.05  % fraction of dissolved P lost with runoff per unit flow [-] % 0.15 or 0.175
beta_1_max = 0.175;
beta_1_vec = linspace(beta_1_min,beta_1_max,num_sets);
beta_2_min = 0.000125; % 000005% fraction of particulate P lost with runoff per unit flow [-] (for 1 m depth, Q_bar=4, 0.00025 -> 1.0 mm loss) % 0.00015 or 0.000175
beta_2_max = 0.000175;
beta_2_vec = linspace(beta_2_min,beta_2_max,num_sets);


% params = lhsdesign(num_sets,2);
u_min = 0; % "natural" P weathering rate % 7 if including subsoil supply here
u_max = 30 ; % max 30 kg/ha fertilizer input % +7 if including subsoil supply here
x1_init = 200; % 226 for 20 cm, 200 for 17.75 cm (7 inches)


%% Run code and identify hitting times for all parameterizations
hit_vector_opt = zeros(num_sets,1);
hit_vector_stquo = zeros(num_sets,1);
hit_vector_constrained = zeros(num_sets,1);
Yr45Reach = repmat(polyshape, 1, num_sets);
horizon = 45; % number of seasons or years to simulate
f = waitbar(0, 'Starting scenarios');
%%%%%debug begin 
% for zz=91:91
%%%%% debug end
for zz=1:num_sets
    % Parameter value, tested for sensitivity
    beta_1 = beta_1_vec(zz);
%     beta_2 = beta_2_vec(zz);
%     beta_1 = beta_1_min + (params(zz,1))*(beta_1_max-beta_1_min);
%     beta_2 = beta_2_min + (params(zz,2))*(beta_2_max-beta_2_min);

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
    
    run('SimpleModel_safe_interv.m'); 
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
%%
toc

%% Plots
figure; plot(beta_1_vec,hit_vector_opt)
xlim([beta_1_min beta_1_max])
title('\beta_1', 'fontsize', 26)
xlabel('Parameter value', 'fontsize', 26)
ylabel('Minimum hitting time (years)', 'fontsize', 26)
ax = gca; % axes handle
ax.FontSize = 18; 
ylim([18 42])
% figure; plot(beta_2_vec,hit_vector_opt)
% xlim([beta_2_min beta_2_max])
% xticks([0.00013 0.00014 0.00015 0.00016 0.00017])
% xticklabels({'0.13','0.14','0.15','0.16','0.17'})
% title('\beta_2', 'fontsize', 26)
% xlabel('Parameter value x 10^3', 'fontsize', 26)
% ylabel('Minimum hitting time (years)', 'fontsize', 26)
% ax = gca; % axes handle
% ax.FontSize = 18; 
% ylim([18 42])

figure; plot(beta_1_vec,hit_vector_stquo)
xlim([beta_1_min beta_1_max])
title('\beta_1', 'fontsize', 26)
xlabel('Parameter value', 'fontsize', 26)
ylabel('Status quo hitting time (years)', 'fontsize', 26)
ax = gca; % axes handle
ax.FontSize = 18; 
ylim([243 423])
% figure; plot(beta_2_vec,hit_vector_stquo)
% xlim([beta_2_min beta_2_max])
% xticks([0.00013 0.00014 0.00015 0.00016 0.00017])
% xticklabels({'0.13','0.14','0.15','0.16','0.17'})
% title('\beta_2', 'fontsize', 26)
% xlabel('Parameter value x 10^3', 'fontsize', 26)
% ylabel('Status quo hitting time (years)', 'fontsize', 26)
% ax = gca; % axes handle
% ax.FontSize = 18; 
% ylim([243 423])