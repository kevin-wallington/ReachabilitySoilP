%% Initialize vectors and set controls
horizon = 750;
X = zeros(horizon, 2); % state: dissolved+adsorbed, internally fixed, streambed
% Initial x1  = 52 g/Mg (Bray P ave) * 1.5 (convert to total ads+dis) * 1.45 Mg/m^3 * 10,000 m^2/ha = 1131 kg/ha (assumes modelling 1 m depth)(226 for .2 m)
X(1,:)=[ x1_init_temp; x2_init_temp]; % estimated IC
O = zeros(horizon, 2); % output: adsorbed pool, river export
hit_time = 0;
switch_time = 0;
singular_time = 0;
fert_traj = zeros(horizon, 1);

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
% pgon=Reach_array(1);
% [r1x,r1y] = boundary(pgon); 
% pgon=Reach_array(2);
% [r2x,r2y] = boundary(pgon); 
% pgon=Reach_array(3);
% [r3x,r3y] = boundary(pgon); 
%% Simulate and plot sample trajectory
fert_rate = u_min; % fertilizer rate [kg/ha] 
for ii=2:horizon
    f_nonlin = [(-alpha_4*(X(ii-1,1)/(X(ii-1,1)+alpha_5))); 0];
    b = [fert_rate;0]; 
    fert_traj(ii) = b(1);
    X(ii,:) = (A * X(ii-1,:)' + f_nonlin + b)'; % update state
    O(ii-1,:) = C * X(ii-1,:)';  
    % singular
    if(inpolygon(X(ii,1),X(ii,2),sx,sy)) % in singular region
        singular_time = ii-1;
        x1line = [X(ii,1)-100 , X(ii,1)+100]; % identify -c3/c4 line
        x2line = [X(ii,2)+(A(2,1)/A(2,2))*(100) , X(ii,2)-(A(2,1)/A(2,2))*(100)];
        if(inpolygon(X(ii,1),X(ii,2),s1x,s1y)) 
%             [x1i,x2i] = polyxpoly(x1line,x2line,s1x,s1y); % find point along line at u-max edge of singular region
            [x1i,x2i] = polyxpoly(x1line,x2line,scon1x,scon1y); % find point along line at u-max edge of singular region - editted 08/09/23 so uses new "singular_control" array, doesn't subtract out prior singular regions
%             x1_bdry = min(x1i); x2_bdry=max(x2i); % guide to right boundary of next target
            x1_bdry = max(x1i); x2_bdry=min(x2i); % guide to left boundary of next target
        end
        if(inpolygon(X(ii,1),X(ii,2),s2x,s2y)) 
%             [x1i,x2i] = polyxpoly(x1line,x2line,s2x,s2y); % find point along line at u-max edge of singular region
               [x1i,x2i] = polyxpoly(x1line,x2line,scon2x,scon2y); % replaces above 08/08/23
%             x1_bdry = min(x1i); x2_bdry=max(x2i); % guide to right boundary of next target
            x1_bdry = max(x1i); x2_bdry=min(x2i); % guide to left boundary of next target
        end
        if(inpolygon(X(ii,1),X(ii,2),s3x,s3y)) 
%             [x1i,x2i] = polyxpoly(x1line,x2line,s3x,s3y); % find point along line at u-max edge of singular region
               [x1i,x2i] = polyxpoly(x1line,x2line,scon3x,scon3y); % replaces above 08/08/23
%             x1_bdry = min(x1i); x2_bdry=max(x2i); % guide to right boundary of next target
            x1_bdry = max(x1i); x2_bdry=min(x2i); % guide to left boundary of next target
        end
        if(inpolygon(X(ii,1),X(ii,2),s4x,s4y)) 
%             [x1i,x2i] = polyxpoly(x1line,x2line,s4x,s4y); % find point along line at u-max edge of singular region
               [x1i,x2i] = polyxpoly(x1line,x2line,scon4x,scon4y); % replaces above 08/08/23
%             x1_bdry = min(x1i); x2_bdry=max(x2i); % guide to right boundary of next target
            x1_bdry = max(x1i); x2_bdry=min(x2i); % guide to left boundary of next target
        end
        if(inpolygon(X(ii,1),X(ii,2),s5x,s5y)) 
%             [x1i,x2i] = polyxpoly(x1line,x2line,s5x,s5y); % find point along line at u-max edge of singular region
               [x1i,x2i] = polyxpoly(x1line,x2line,scon5x,scon5y); % replaces above 08/08/23
%             x1_bdry = min(x1i); x2_bdry=max(x2i); % guide to right boundary of next target
            x1_bdry = max(x1i); x2_bdry=min(x2i); % guide to left boundary of next target
        end        
        if(inpolygon(X(ii,1),X(ii,2),s6x,s6y)) 
%             [x1i,x2i] = polyxpoly(x1line,x2line,s6x,s6y); % find point along line at u-max edge of singular region
               [x1i,x2i] = polyxpoly(x1line,x2line,scon6x,scon6y); % replaces above 08/08/23
%             x1_bdry = min(x1i); x2_bdry=max(x2i); % guide to right boundary of next target
            x1_bdry = max(x1i); x2_bdry=min(x2i); % guide to left boundary of next target
                end   
        if(inpolygon(X(ii,1),X(ii,2),s7x,s7y)) 
%             [x1i,x2i] = polyxpoly(x1line,x2line,s7x,s7y); % find point along line at u-max edge of singular region
               [x1i,x2i] = polyxpoly(x1line,x2line,scon7x,scon7y); % replaces above 08/08/23
%             x1_bdry = min(x1i); x2_bdry=max(x2i); % guide to right boundary of next target
            x1_bdry = max(x1i); x2_bdry=min(x2i); % guide to left boundary of next target
        end   
        u_diff = A(1,1)*(x1_bdry-X(ii,1)) + A(1,2)*(x2_bdry-X(ii,2)) - ...
           alpha_4*(x1_bdry/(x1_bdry+alpha_5) - X(ii,1)/(X(ii,1)+alpha_5)); % calculate delta-u
%         fert_rate = u_max + u_diff; % for guiding to right boundary, u_diff will be negative
        fert_rate = u_min + u_diff; % for guiding to left boundary, u_diff will be positive
        if (fert_rate <= u_min); fert_rate = u_min; end
        if (fert_rate >= u_max); fert_rate = u_max; end
        b = [fert_rate;0];
    end
    if(inpolygon(X(ii,1),X(ii,2), ux,uy) && hit_time == 0) % in region to switch to fert max
        if switch_time == 0
            switch_time = ii-1;
        end
        fert_rate = u_max;
        b = [fert_rate;0];
    end
    if(inpolygon(X(ii,1),X(ii,2), cx,cy) && hit_time == 0) % in target region for first time
        hit_time = ii-1;
    end
    if(inpolygon(X(ii,1),X(ii,2), cx,cy)) % in target region
        f_nonlin = [(-alpha_4*(X(ii,1)/(X(ii,1)+alpha_5))); 0];
        fert_rate = 17 + u_min; % draw down a bit
        b = [fert_rate;0];
        X1_temp = (A(1,:) * X(ii,:)' + f_nonlin(1) + b(1))'; % next state estimate
%         if(X1_temp <= (X(ii,1)))
       fert_rate= fert_rate + X(ii,1) - X1_temp; 
       if (fert_rate <= u_min); fert_rate = u_min; end
       if (fert_rate >= u_max); fert_rate = u_max; end
       b = [fert_rate; 0]; 
%         end
    end
end

%% Plot state space evolution
% figure
% hold on
% plot(SafeInv,'FaceColor',[0.5 0.5 0.5], 'FaceAlpha',1); 
% plot(Singular_array_all,'FaceColor',[0.98 0.85 0.54], 'FaceAlpha',1); 
% plot(X(1:ii,1), X(1:ii,2),'-', 'color', "#77AC30", 'linewidth', 3) %"#77AC30" for optimal trajectory
% plot(X(1:horizon,1), X(1:horizon,2),'-', 'color', "#77AC30", 'linewidth', 3) %"#77AC30" for optimal trajectory
% % plot formatting
% xlim([0 500])
% ylim([000 10000])
% ax = gca; % axes handle
% ax.FontSize = 18; 
% title('', 'fontsize', 26)
% xlabel('Short-term P (kg/ha)', 'fontsize', 26)
% ylabel('Long-term P (kg/ha)', 'fontsize', 26)
% % legend('', 'Location', 'Northeast')
% % legend boxoff
% % plot initial point
% plot([x1_init],[x2_init],'pentagram','MarkerFaceColor','black','MarkerEdgeColor','black', 'MarkerSize',20,'linewidth',2);

