% Variable interim constraint on nutrient deficit
Def_lim = 0.75; % nutrient deficit limit (i.e. 75% of max)
x1_lim = Def_lim/(1-Def_lim)*alpha_5; % overwrite below for exploring scenarios
scenarios = 56; % was 61
x1_lim = linspace(50,105,scenarios); % was 60 to 120

%% Initialize vectors and set controls
horizon = 100;
X = cell(scenarios,1);
O = cell(scenarios,1);
X(:) = {zeros(horizon, 2)}; % state: dissolved+adsorbed, internally fixed, streambed
% Initial x1  = 52 g/Mg (Bray P ave) * 1.5 (convert to total ads+dis) * 1.45 Mg/m^3 * 10,000 m^2/ha = 1131 kg/ha (assumes modelling 1 m depth)(226 for .2 m)
for ii=1:scenarios%+1
    X{ii}(1,:)=[ x1_init; x2_init]; % estimated IC
end
O(:) = {zeros(horizon, 2)}; % output: adsorbed pool, river export
hit_time = zeros(scenarios,1);
switch_time = zeros(scenarios,1);
singular_time = zeros(scenarios,1);
fert_traj = cell(scenarios,1);
fert_traj(:) = {zeros(horizon, 1)};

%% Target and switching region
pgon=SafeInv; % Safe_Inv created by SimpleModel_Invariance.m
polyout = polybuffer(pgon,0.1);
[cx,cy] = boundary(polyout); 
pgon=umax_region_pgon;
polyout = polybuffer(pgon,0.01);
[ux,uy] = boundary(polyout); 
pgon=Singular_array_all;
[sx,sy] = boundary(pgon); 
pgon=Singular_array(1);
[s1x,s1y] = boundary(pgon); 
pgon=Singular_array(2);
[s2x,s2y] = boundary(pgon); 
pgon=Singular_array(3);
[s3x,s3y] = boundary(pgon); 
pgon=Singular_array(4);
[s4x,s4y] = boundary(pgon);
pgon=Singular_array(5);
[s5x,s5y] = boundary(pgon);
pgon=Singular_array(6);
[s6x,s6y] = boundary(pgon);
pgon=Singular_array(7);
[s7x,s7y] = boundary(pgon);
pgon=Reach_array(1);
[r1x,r1y] = boundary(pgon); 
pgon=Reach_array(2);
[r2x,r2y] = boundary(pgon); 
pgon=Reach_array(3);
[r3x,r3y] = boundary(pgon); 

pgon=Singular_array_control(1);
[scon1x,scon1y] = boundary(pgon); 
pgon=Singular_array_control(2);
[scon2x,scon2y] = boundary(pgon); 
pgon=Singular_array_control(3);
[scon3x,scon3y] = boundary(pgon); 
pgon=Singular_array_control(4);
[scon4x,scon4y] = boundary(pgon);
pgon=Singular_array_control(5);
[scon5x,scon5y] = boundary(pgon);
pgon=Singular_array_control(6);
[scon6x,scon6y] = boundary(pgon);
pgon=Singular_array_control(7);
[scon7x,scon7y] = boundary(pgon);

%% Simulate and plot sample trajectory
for kk=1:scenarios
    fert_rate = u_min; % fertilizer rate [kg/ha] 
    for ii=2:horizon
        f_nonlin = [(-alpha_4*(X{kk}(ii-1,1)/(X{kk}(ii-1,1)+alpha_5))); 0];
        b = [fert_rate;0]; 
        X1_temp = (A(1,:) * X{kk}(ii-1,:)' + f_nonlin(1) + b(1))'; % update state
        if((X1_temp<=x1_lim(kk)) )%&& singular_time(kk)==0)
           adj_fert_rate=fert_rate+x1_lim(kk)-X1_temp; 
           b = [adj_fert_rate; 0]; 
        end
        fert_traj{kk}(ii) = b(1);
        X{kk}(ii,:) = (A * X{kk}(ii-1,:)' + f_nonlin + b)'; % update state
        O{kk}(ii-1,:) = C * X{kk}(ii-1,:)';  
        % singular
        if(inpolygon(X{kk}(ii,1),X{kk}(ii,2),sx,sy)) % in singular region
            singular_time(kk) = ii-1;
            x1line = [X{kk}(ii,1)-100 , X{kk}(ii,1)+100]; % identify -c3/c4 line
            x2line = [X{kk}(ii,2)+(A(2,1)/A(2,2))*(100) , X{kk}(ii,2)-(A(2,1)/A(2,2))*(100)];
            % guide to left boundary of next target
%             if(inpolygon(X{kk}(ii,1),X{kk}(ii,2),s1x,s1y)) 
%                 [x1i,x2i] = polyxpoly(x1line,x2line,s1x,s1y); % find point along line at u-min edge of singular region
%                 x1_bdry = max(x1i); x2_bdry=min(x2i);
%             end
%             if(inpolygon(X{kk}(ii,1),X{kk}(ii,2),s2x,s2y)) 
%                 [x1i,x2i] = polyxpoly(x1line,x2line,s2x,s2y); % find point along line at u-min edge of singular region
%                 x1_bdry = max(x1i); x2_bdry=min(x2i);
%             end
%             if(inpolygon(X{kk}(ii,1),X{kk}(ii,2),s3x,s3y)) 
%                 [x1i,x2i] = polyxpoly(x1line,x2line,s3x,s3y); % find point along line at u-min edge of singular region
%                 x1_bdry = max(x1i); x2_bdry=min(x2i);
%             end         
            % guide to right boundary of next target
            if(inpolygon(X{kk}(ii,1),X{kk}(ii,2),s1x,s1y)) 
%                 [x1i,x2i] = polyxpoly(x1line,x2line,s1x,s1y); % find point along line at u-max edge of singular region
%                 x1_bdry = min(x1i); x2_bdry=max(x2i);
                [x1i,x2i] = polyxpoly(x1line,x2line,scon1x,scon1y);
                x1_bdry = max(x1i); x2_bdry=min(x2i); % guide to left boundary of next target
            end
            if(inpolygon(X{kk}(ii,1),X{kk}(ii,2),s2x,s2y)) 
%                 [x1i,x2i] = polyxpoly(x1line,x2line,s2x,s2y); % find point along line at u-max edge of singular region
%                 x1_bdry = min(x1i); x2_bdry=max(x2i);
                [x1i,x2i] = polyxpoly(x1line,x2line,scon2x,scon2y);
                x1_bdry = max(x1i); x2_bdry=min(x2i); % guide to left boundary of next target
            end
            if(inpolygon(X{kk}(ii,1),X{kk}(ii,2),s3x,s3y)) 
%                 [x1i,x2i] = polyxpoly(x1line,x2line,s3x,s3y); % find point along line at u-max edge of singular region
%                 x1_bdry = min(x1i); x2_bdry=max(x2i);
                [x1i,x2i] = polyxpoly(x1line,x2line,scon3x,scon3y);
                x1_bdry = max(x1i); x2_bdry=min(x2i); % guide to left boundary of next target
            end
            if(inpolygon(X{kk}(ii,1),X{kk}(ii,2),s4x,s4y)) 
%                 [x1i,x2i] = polyxpoly(x1line,x2line,s4x,s4y); % find point along line at u-max edge of singular region
%                 x1_bdry = min(x1i); x2_bdry=max(x2i);
                [x1i,x2i] = polyxpoly(x1line,x2line,scon4x,scon4y);
                x1_bdry = max(x1i); x2_bdry=min(x2i); % guide to left boundary of next target
            end
            if(inpolygon(X{kk}(ii,1),X{kk}(ii,2),s5x,s5y)) 
%                 [x1i,x2i] = polyxpoly(x1line,x2line,s5x,s5y); % find point along line at u-max edge of singular region
%                 x1_bdry = min(x1i); x2_bdry=max(x2i);
                [x1i,x2i] = polyxpoly(x1line,x2line,scon5x,scon5y);
                x1_bdry = max(x1i); x2_bdry=min(x2i); % guide to left boundary of next target
            end        
            if(inpolygon(X{kk}(ii,1),X{kk}(ii,2),s6x,s6y)) 
%                 [x1i,x2i] = polyxpoly(x1line,x2line,s6x,s6y); % find point along line at u-max edge of singular region
%                 x1_bdry = min(x1i); x2_bdry=max(x2i);
                [x1i,x2i] = polyxpoly(x1line,x2line,scon6x,scon6y);
                x1_bdry = max(x1i); x2_bdry=min(x2i); % guide to left boundary of next target
            end   
            if(inpolygon(X{kk}(ii,1),X{kk}(ii,2),s7x,s7y)) 
%                 [x1i,x2i] = polyxpoly(x1line,x2line,s7x,s7y); % find point along line at u-max edge of singular region
%                 x1_bdry = min(x1i); x2_bdry=max(x2i);
                [x1i,x2i] = polyxpoly(x1line,x2line,scon7x,scon7y);
                x1_bdry = max(x1i); x2_bdry=min(x2i); % guide to left boundary of next target
            end   

            u_diff = A(1,1)*(x1_bdry-X{kk}(ii,1)) + A(1,2)*(x2_bdry-X{kk}(ii,2)) - ...
               alpha_4*(x1_bdry/(x1_bdry+alpha_5) - X{kk}(ii,1)/(X{kk}(ii,1)+alpha_5)); % calculate delta-u
%             fert_rate = u_max + u_diff; % for guiding to right boundary, u_diff will be negative
            fert_rate = u_diff + u_min; % for guiding to left boundary, u_diff will be positive
            if (fert_rate <= u_min); fert_rate = u_min; end
            if (fert_rate >= u_max); fert_rate = u_max; end
            b = [fert_rate;0];
        end
        if(inpolygon(X{kk}(ii,1),X{kk}(ii,2), ux,uy) && hit_time(kk) == 0) % in region to switch to fert max
            if switch_time(kk) == 0
                switch_time(kk) = ii-1;
            end
            fert_rate = 30; % +7 if including subsoil supply here
            b = [fert_rate;0];
        end
        if(inpolygon(X{kk}(ii,1),X{kk}(ii,2), cx,cy) && hit_time(kk) == 0) % in target region for first time
            hit_time(kk) = ii-1;
        end
        if(inpolygon(X{kk}(ii,1),X{kk}(ii,2), cx,cy)) % in target region
            f_nonlin = [(-alpha_4*(X{kk}(ii,1)/(X{kk}(ii,1)+alpha_5))); 0];
            fert_rate = 17; % draw down a bit % +7 if including subsoil supply rate
            b = [fert_rate;0];
            X1_temp = (A(1,:) * X{kk}(ii,:)' + f_nonlin(1) + b(1))'; % next state estimate
%             if(X1_temp <= (X{kk}(ii,1)))
               fert_rate= fert_rate + X{kk}(ii,1) - X1_temp; 
               b = [fert_rate; 0]; 
%             end
        end
    end
end
opt_hit_time = hit_time(1);
constr_hit_time = hit_time(scenarios);
%% status quo trajectory
stquo_traj = zeros(horizon, 2);
stquo_traj(1,:)=[ x1_init; x2_init];
fert_rate = 21; % fertilizer rate [kg/ha] 
horizon_stquo = 500;
stquo_hit_time = 0;
for ii=2:horizon_stquo
    f_nonlin = [(-alpha_4*(stquo_traj(ii-1,1)/(stquo_traj(ii-1,1)+alpha_5))); 0];
    b = [fert_rate;0]; 
    stquo_traj(ii,:) = (A * stquo_traj(ii-1,:)' + f_nonlin + b)'; % update state 
    if(inpolygon(stquo_traj(ii,1),stquo_traj(ii,2), cx,cy) && stquo_hit_time == 0) % in target region for first time
        stquo_hit_time = ii-1;
    end
end

%% plot target polygon for deterministic x1,x2, initial point, equilibrium b/w x1-x2
figure
hold on

x1_min = Def_lim/(1-Def_lim)*alpha_5;
x2_max = (Load_lim - Q_bar*(beta_1*alpha_2+beta_2*(1-alpha_2))*x1_min) / (Q_bar*beta_2);
x1_max = Load_lim /(Q_bar*(beta_1*alpha_2+beta_2*(1-alpha_2)));
A = [x1_min, x1_max];
B = [x2_max 0];
xlimits = [0 220];
m = (B(2)-B(1))/(A(2)-A(1));
n = -A(2)*m + B(2);
y1 = m*xlimits(1) + n;
y2 = m*xlimits(2) + n;

x2 = [xlimits, fliplr(xlimits)];
inBetween = [[10000; 10000]', fliplr([y1 y2])];
fill(x2, inBetween, [0.93,0.69,0.13], 'FaceAlpha', 0.1, 'EdgeColor', 'none');

% And a line for the crop nutrient deficit
xlimits = [0 x1_min];
x2 = [xlimits, fliplr(xlimits)];
highline = 200*ones(scenarios,1);
inBetween = [[10000; 10000]', fliplr([0 0])];
fill(x2, inBetween, [0.42,0.12,0.12], 'FaceAlpha', 0.1, 'EdgeColor', 'none');

% plot(pmask,'FaceColor',[0 0.2 0.35],'FaceAlpha',1);
% plot(umax_region_pgon,'FaceColor',[0.42 0.58 0.69],'FaceAlpha',1);
plot(SafeSet, 'FaceColor',[0 0 0], 'FaceAlpha',1);
plot(SafeInv,'FaceColor',[0.5 0.5 0.5], 'FaceAlpha',1); 
% plot(Singular_array_all,'FaceColor',[0.98 0.85 0.54], 'FaceAlpha',1); 
% lgd = legend('Minimum control region','Maximum control region','Target','Singular region', 'Location', 'Northwest', 'fontsize', 26);
% lgd.Color = 'white';
% legend boxon
% lgd.AutoUpdate = 'off';
%% Plot state space evolution
plot(X{scenarios}(1:constr_hit_time+1,1), X{scenarios}(1:constr_hit_time+1,2),'-', 'color', "black", 'linewidth', 6) %constrained trajectory
plot(X{1}(1:opt_hit_time+1,1), X{1}(1:opt_hit_time+1,2),'-', 'color', "#77AC30", 'linewidth', 6) %"#77AC30" for optimal trajectory
% for kk=[1,scenarios]
%     plot(X{kk}(1:horizon,1), X{kk}(1:horizon,2),'-', 'color', "#77AC30", 'linewidth', 4) %"#77AC30" for optimal trajectory
% end
plot(stquo_traj(1:stquo_hit_time+1,1), stquo_traj(1:stquo_hit_time+1,2),'-', 'color', "black", 'linewidth', 6, 'LineStyle', ':') %status quo trajectory
% plot([X{1}(31,1)], [X{1}(31,2)],'*','MarkerFaceColor','white','MarkerEdgeColor',[0.75 0 0.15], 'MarkerSize',25,'linewidth',6); % old-27, new-??
% plot([X{1}(40,1)], [X{1}(40,2)],'o','MarkerFaceColor','black','MarkerEdgeColor','black', 'MarkerSize',10,'linewidth',6); % old-39, new-43
% plot([X{1}(41,1)], [X{1}(41,2)],'o','MarkerFaceColor','black','MarkerEdgeColor','black', 'MarkerSize',10,'linewidth',6); % old-40, new-44
% plot([X{1}(42,1)], [X{1}(42,2)],'o','MarkerFaceColor','black','MarkerEdgeColor','black', 'MarkerSize',10,'linewidth',6); % old-41, new-none
% plot formatting
xlim([30 220])
ylim([mask_y 1800])
ax = gca; % axes handle
ax.FontSize = 26; 
title('', 'fontsize', 26)
xlabel('Labile phosphorus (kg/ha)', 'fontsize', 36)
ylabel('Stable phosphorus (kg/ha)', 'fontsize', 36)
% legend('', 'Location', 'Northwest')
% legend boxoff
% plot initial point
% plot([x1_init],[x2_init],'pentagram','MarkerFaceColor','black','MarkerEdgeColor','black', 'MarkerSize',20,'linewidth',2);
plot([x1_init],[x2_init],'pentagram','MarkerFaceColor',[0.75 0 0.15],'MarkerEdgeColor','white', 'MarkerSize',100,'linewidth',2);

%% Plot hitting time
figure 
set(gcf,'position',[10,10,1650,1200])
hold on
deficit = 105 - x1_lim;
x2 = [deficit(1:scenarios), fliplr(deficit(1:scenarios))];
% x2 = [x1_lim(1:scenarios), fliplr(x1_lim(1:scenarios))];
highline = 100*ones(scenarios,1);
inBetween = [highline(1:scenarios)', fliplr(hit_time(1:scenarios)')];
fill(x2, inBetween, [0.94,0.94,0.94]);

grad1 = zeros(17,3);
for k = 1:17, grad1(k,:) = [0.34,0.34,0.99] ; end
grad2 = colorGradient([0 0 0],[0.34,0.34,0.99],36);
grad = [grad2' grad1']';
% for k = 1:17, grad1(k,:) = [0.34,0.34,0.99] ; end
% grad2 = colorGradient([0.34,0.34,0.99],[0 0 0],36);
% grad2 = colorGradient([0.47,0.67,0.19],[0 0 0],36);
% grad = [grad1' grad2']';
xdata = deficit(1:scenarios);
% xdata = x1_lim(1:scenarios);
ydata = hit_time(1:scenarios)';
zdata = zeros(size(X))';
% zdata = zeros(size(X));
col = xdata;  % This is the color, vary with x in this case.
surface([xdata;xdata],[ydata;ydata],[zdata;zdata],[col;col],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',6);
colormap(grad);

% plot(x1_lim(1:scenarios), hit_time(1:scenarios),'linewidth',6, 'Color', [0.47,0.67,0.19]);
xlim([0 55])
% xlim([50 105])
ylim([42 66])
ax = gca; % axes handle
ax.FontSize = 26; 
title('', 'fontsize', 26)
xlabel('Interim crop nutrient deficit [kg/ha]', 'fontsize', 36)
% xlabel('Minimum interim labile P level [kg/ha]', 'fontsize', 36)
ylabel('Time to joint target (years)', 'fontsize', 36)
plot([1], [43],'pentagram','MarkerFaceColor',[0.34,0.34,0.99],'MarkerEdgeColor','black', 'MarkerSize',75,'linewidth',6); % old-27, new-??
% plot([104], [43],'pentagram','MarkerFaceColor',[0.34,0.34,0.99],'MarkerEdgeColor','black', 'MarkerSize',75,'linewidth',6); % old-27, new-??



