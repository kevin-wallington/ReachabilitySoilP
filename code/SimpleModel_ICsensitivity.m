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
x1_init_vec = (150:5:250); % 21 ICs
x2_diseq_vec = (-410.5:41.05:410.5); % 21 ICs
num_sets = length(x1_init_vec)*length(x2_diseq_vec);
u_min = 0;%7; % "natural" P weathering rate
u_max = 30;% + 7; % max 30 kg/ha fertilizer input

%% Run code and identify hitting times for all parameterizations
hit_vector_opt = zeros(num_sets,1);
hit_vector_stquo = zeros(num_sets,1);
hit_vector_constrained = zeros(num_sets,1);
hitmap = zeros(length(x2_diseq_vec),length(x1_init_vec));
% x1_plot = zeros(num_sets,1);
% x2_plot = zeros(num_sets,1);
Yr45Reach = repmat(polyshape, 1, num_sets);
horizon = 45; % number of seasons or years to simulate
f = waitbar(0, 'Starting scenarios');
%%%%%debug begin 
% for zz=305:305
%%%%% debug end
for zz=1:num_sets
    x1_id = (floor((zz-1)/length(x2_diseq_vec))+1);
    x1_init = x1_init_vec(x1_id); 
    x2_id =  mod(zz,length(x2_diseq_vec));
    if(x2_id==0); x2_id= length(x2_diseq_vec); end
    x2_init = x1_init/alpha_1 + x2_diseq_vec(x2_id);
    Q_bar = Loss/((beta_1*alpha_2+beta_2*(1-alpha_2))*x1_init + beta_2*x2_init);
%     x1_plot(zz) = x1_init;
%     x2_plot(zz)= x2_init;

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
        hitmap(length(x2_diseq_vec)+1-x2_id,x1_id) = hit_time;
        run('SimpleModel_stquotime.m'); % calculate hitting time for no intermediate contstraint
        hit_vector_stquo(zz) = hit_time;
        run('SimpleModel_constrainedtime.m'); % calculate hitting time for no intermediate contstraint
        hit_vector_constrained(zz) = hit_time;
    else
        hit_vector_opt(zz) = NaN;
        hitmap((length(x2_diseq_vec)+1-x2_id),x1_id) = NaN;
        hit_vector_stquo(zz) = NaN;
        hit_vector_constrained(zz) = NaN;
    end 
    waitbar(zz/num_sets, f, sprintf('Progress: %d %%', floor(zz/num_sets*100)));
end
close(f)
%%
toc

%% Plots
% figure; hist(hit_vector_opt(find(hit_vector_opt>0)),9)
% figure; hist(hit_vector_constrained(find(hit_vector_constrained>0)),40)
% figure; hist(hit_vector_stquo(find(hit_vector_stquo>0)),40)

figure
h=heatmap(hitmap,'CellLabelColor','none','Colormap',parula,'ColorLimits',[min(hit_vector_opt) max(hit_vector_opt)],...
    'MissingDataColor',[0.9 0.9 0.9],'GridVisible','off');
xlabels = num2str(x1_init_vec');
ylabels = num2str(flip(x2_diseq_vec'));
h.XDisplayLabels = xlabels;
h.YDisplayLabels = ylabels;
% scatter(x1_plot, x2_plot, 25,hit_vector_opt, 'filled');
% xlim([150 250])
% ylim([1100 2100])
% ax = gca; % axes handle
% ax.FontSize = 18; 
% title('', 'fontsize', 26)
% xlabel('Short-term P (kg/ha)', 'fontsize', 26)
% ylabel('Long-term P (kg/ha)', 'fontsize', 26)
% a = colorbar;
% a.Label.String = 'Minimum hitting time';
% plot([200],[1642],'pentagram','MarkerFaceColor','black','MarkerEdgeColor','black', 'MarkerSize',20,'linewidth',2);