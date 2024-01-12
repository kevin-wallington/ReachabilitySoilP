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
Loss = 2.0; % [kg/ha] average annual losses currently
x1_init = 200; 
x2_init = x1_init/alpha_1 + 42;
umin_vec = (3:1:10); % 8 values
fertmax_vec = (25:1:35); % 11 values
num_sets = length(umin_vec)*length(fertmax_vec);
% u_min = 7; % "natural" P weathering rate
% u_max = 30 + 7; % max 30 kg/ha fertilizer input

%% Run code and identify hitting times for all parameterizations
hit_vector_opt = zeros(num_sets,1);
hit_vector_stquo = zeros(num_sets,1);
hit_vector_constrained = zeros(num_sets,1);
hitmap = zeros(length(fertmax_vec),length(umin_vec));
hitmap_constr = zeros(length(fertmax_vec),length(umin_vec));
Yr45Reach = repmat(polyshape, 1, num_sets);
horizon = 45; % number of seasons or years to simulate
f = waitbar(0, 'Starting scenarios');
%%%%%debug begin 
% for zz=41:41
%%%%% debug end
for zz=1:num_sets
    min_id = (floor((zz-1)/length(fertmax_vec))+1);
    u_min = umin_vec(min_id); 
    max_id =  mod(zz,length(fertmax_vec));
    if(max_id==0); max_id= length(fertmax_vec); end
    u_max = u_min + fertmax_vec(max_id);

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
        hit_vector_opt(zz) = hit_time;
        hitmap((length(fertmax_vec)+1-max_id),min_id) = hit_time;
        run('SimpleModel_stquotime.m'); % calculate hitting time for no intermediate contstraint
        hit_vector_stquo(zz) = hit_time;
        run('SimpleModel_constrainedtime.m'); % calculate hitting time for no intermediate contstraint
        hit_vector_constrained(zz) = hit_time;
        hitmap_constr((length(fertmax_vec)+1-max_id),min_id) = hit_time;
    else
        hit_vector_opt(zz) = NaN;
        hitmap((length(fertmax_vec)+1-max_id),min_id) = Nan;
        hit_vector_stquo(zz) = NaN;
        hit_vector_constrained(zz) = NaN;
    end 
    waitbar(zz/num_sets, f, sprintf('Progress: %d %%', floor(zz/num_sets*100)));
end
close(f)
%%
toc

%% Plots
% figure; hist(hit_vector_opt(find(hit_vector_opt>0)),10)
% figure; hist(hit_vector_constrained(find(hit_vector_constrained>0)),40)
% figure; hist(hit_vector_stquo(find(hit_vector_stquo>0)),40)

figure
h=heatmap(hitmap,'CellLabelColor','none','Colormap',parula,'ColorLimits',[min(hit_vector_opt) max(hit_vector_opt)],...
    'MissingDataColor',[0.9 0.9 0.9],'GridVisible','off');
xlabels = num2str(umin_vec');
ylabels = num2str(flip(fertmax_vec'));
h.XDisplayLabels = xlabels;
h.YDisplayLabels = ylabels;

% figure
% h=heatmap(hitmap_constr,'CellLabelColor','none','Colormap',parula,'ColorLimits',[min(hit_vector_constrained) max(hit_vector_constrained)],...
%     'MissingDataColor',[0.9 0.9 0.9],'GridVisible','off');
% xlabels = num2str(umin_vec');
% ylabels = num2str(flip(fertmax_vec'));
% h.XDisplayLabels = xlabels;
% h.YDisplayLabels = ylabels;