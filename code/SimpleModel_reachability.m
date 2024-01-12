%% Create target polygon points vector and mask for state space
pgon= SafeInv;
umaxback_pgon = SafeInv;
[cx,cy] = boundary(SafeInv);
mask_y= min(cy); % DEBUG 07/20/23 - cuts off bottom of polygons so don't have issues
mask_x = 30;
% mask_y = 0;
mask_points = [mask_x mask_y;mask_x 2500;1500 2500;1500 mask_y];
pmask = polyshape(mask_points);

%% Back propagate
horizon = 50; % number of seasons or years to simulate
Singular_array = repmat(polyshape, 1, horizon);
Singular_array_control = repmat(polyshape, 1, horizon);
BBreach_array = repmat(polyshape, 1, horizon);
Reach_array = repmat(polyshape, 1, horizon);
for jj=1:horizon % back propagation loop
    % sample boundary of polygon from last step
    samples = 0;
    other_poly_points=[];
    [cx,cy] = boundary(pgon); % will wrap back to original pt to close polygon
    num_pts = length(cx);
    for kk=1:(num_pts-1) %generate points along boundary of polygon
       d=pdist([cx(kk) cy(kk); cx(kk+1) cy(kk+1)]);
       seg_pts = floor(d/10)+2;
       x_seg = linspace(cx(kk),cx(kk+1),seg_pts);
       y_seg = linspace(cy(kk),cy(kk+1),seg_pts);
       x_seg(end)=[]; y_seg(end)=[];
       new_pts = [x_seg; y_seg];
       other_poly_points = [other_poly_points' new_pts]';
    end    
    other_poly_points(abs(other_poly_points)<0.1)=0; % otherwise, points very close to zero
    % For u_max only, b/c can't go max->min or max->singular
    maxback_poly_points=[];
    [mbx,mby] = boundary(umaxback_pgon); % will wrap back to original pt to close polygon
    num_pts = length(mbx);
    for kk=1:(num_pts-1) %generate points along boundary of polygon
       d=pdist([mbx(kk) mby(kk); mbx(kk+1) mby(kk+1)]);
       seg_pts = floor(d/10)+2;
       x_seg = linspace(mbx(kk),mbx(kk+1),seg_pts);
       y_seg = linspace(mby(kk),mby(kk+1),seg_pts);
       x_seg(end)=[]; y_seg(end)=[];
       new_pts = [x_seg; y_seg];
       maxback_poly_points = [maxback_poly_points' new_pts]';
    end    
    maxback_poly_points(abs(maxback_poly_points)<0.1)=0; % otherwise, points very close to zero
    % Backpropogate with bang-bang control
    umax_poly_points = zeros(samples,2);
    umin_poly_points = zeros(samples,2);
    for ll=1:2
        % eliminate = [];
        poly_points=[];
        if(ll==1)
            u=u_max; 
            poly_points=maxback_poly_points; 
        end
        if(ll==2)
            u=u_min; 
            poly_points=other_poly_points; 
        end
        samples=length(poly_points);
        for ii=1:samples % back propagate each point
            c5 = A(1,1) - A(1,2) * A(2,1) / A(2,2);
            c6 = (A(1,2) / A(2,2)) * poly_points(ii,2) + u;
            a = c5; 
            b = c5*alpha_5 - alpha_4 + c6 - poly_points(ii,1);
            c = c6*alpha_5 - alpha_5*poly_points(ii,1);
            root_out = roots([a b c]);
            if (~isreal(root_out))
                % eliminate = [eliminate, ii]; 
                pup = alpha_4*poly_points(ii,1)/(poly_points(ii,1)+alpha_5);
                d_pup = alpha_4*alpha_5/(poly_points(ii,1)+alpha_5)^2;
                if(d_pup>=0.95*A(1,1));d_pup=0.95*A(1,1);end % to prevent instability near x1 = 0
                b1 = u - pup + d_pup*poly_points(ii,1);
                b = [b1;0]; 
                back_prop = poly_points(ii,:)' - b;
                A_loc = A;
                A_loc(1,1)=A_loc(1,1)-d_pup;
                init_est = (A_loc\back_prop)';
                new_x1 = init_est(1);
            else
                new_x1 = max(root_out);
                % max_minx1 = max(max_minx1, min(root_out));
            end
            new_x2 = poly_points(ii,2) / A(2,2) -  (A(2,1)/A(2,2))*new_x1;
            if(ll==1); umax_poly_points(ii,:) = [new_x1 new_x2]; end
            if(ll==2); umin_poly_points(ii,:) = [new_x1 new_x2]; end 
        end
        % if(ll==1); umax_poly_points(eliminate,:) = []; end
        % if(ll==2); umin_poly_points(eliminate,:) = []; end
    end
    umax_poly_points(abs(umax_poly_points)<0.1)=0;
    umin_poly_points(abs(umin_poly_points)<0.1)=0;
    umax_pgon=polyshape(umax_poly_points);
    umin_pgon=polyshape(umin_poly_points);
    new_pgon = union(umax_pgon, umin_pgon); % NEW
%     new_pgon = union(pgon,umax_pgon);
%     new_pgon = union(new_pgon,umin_pgon);
    new_pgon = intersect(new_pgon,pmask);
    umaxback_pgon = union(umax_pgon,SafeInv); %%%%%%%%%%%%%%     replace all umaxback_pgon with umax_pgon to revert to earlier version
    umaxback_pgon = intersect(umaxback_pgon,pmask); %%%%%%%%%%%%%%
    BBreach_array(jj)=new_pgon;
    [Rch_cx, Rch_cy] = convU(BBreach_array(jj), A, mask_x, mask_y);  
%     pgon = polyshape(Rch_cx, Rch_cy); 
    new_pgon = polyshape(Rch_cx, Rch_cy); % NEW
    pgon = union(pgon,new_pgon); %NEW
    Reach_array(jj)=pgon;
    Singular_array(jj) = subtract(Reach_array(jj),BBreach_array(jj));
    Singular_array(jj) = subtract(Singular_array(jj),SafeInv);
    Singular_array_control(jj) = Singular_array(jj);
    if(jj>=2);Singular_array(jj) = subtract(Singular_array(jj),Singular_array(jj-1));end
end

%% Singular region excludes bang-bang reachable region
Singular_array_all = Singular_array(1);
for ii=2:4 % small regions due to numerical errors at 7 steps (4 is enough for meaningful regions in standard params)
    Singular_array_all = union(Singular_array_all,Singular_array(ii));
end
Singular_array_all = subtract(Singular_array_all,SafeInv);

%% Umaxback_pgon should not include regions to right of target
% x1_min = Def_lim/(1-Def_lim)*alpha_5;
% x2_max = (Load_lim - Q_bar*(beta_1*alpha_2+beta_2*(1-alpha_2))*x1_min) / (Q_bar*beta_2);
% x1_max = Load_lim /(Q_bar*(beta_1*alpha_2+beta_2*(1-alpha_2)));
% cx = [x1_min, x1_min, 10000, x1_max]; % corner points, x-coord
% cy = [x2_max, 5000, 0, 0]; % corner points, y-coord
% poly_points = [cx; cy]';
% ignore_pgon=polyshape(poly_points);
% umaxback_plot = subtract(umaxback_pgon,ignore_pgon);
% background_points = [mask_x mask_y;mask_x 600;100 600;100 mask_y];
% pgon = polyshape(background_points);
% umaxback_plot = union(umaxback_plot, pgon);

%% Plots
% Reachable sets
figure
hold on
plot(Reach_array(50),'FaceColor',[0.8 0.8 1.0],'EdgeColor',[1 1 1],'FaceAlpha',1);
plot(Reach_array(42),'FaceColor',[0.7 0.7 0.9],'EdgeColor',[1 1 1],'FaceAlpha',1);
plot(Reach_array(10),'FaceColor',[0.4 0.4 0.7],'EdgeColor',[1 1 1],'FaceAlpha',1); % 0 0.2 0.35
plot(Reach_array(3),'FaceColor',[0.2 0.2 0.5],'EdgeColor',[1 1 1],'FaceAlpha',1); % [0.42 0.58 0.69]
plot(Reach_array(2),'FaceColor',[0.1 0.1 0.3],'EdgeColor',[1 1 1],'FaceAlpha',1); % [0 0.2 0.35]
plot(Reach_array(1),'FaceColor',[0 0.0 0.1],'EdgeColor',[1 1 1],'FaceAlpha',1); % [0.42 0.58 0.69]
plot(SafeInv,'FaceColor',[0.5 0.5 0.5],'EdgeColor',[1 1 1],'FaceAlpha',1);
% plot initial condition estimate
plot([x1_init],[x2_init],'pentagram','MarkerFaceColor',[0.75 0 0.15],'MarkerEdgeColor','white', 'MarkerSize',50,'linewidth',2);
% plot formatting
xlim([mask_x 220])
ylim([mask_y 1800])
ax = gca; % axes handle
ax.FontSize = 18; 
title('', 'fontsize', 26)
xlabel('Labile P (kg/ha)', 'fontsize', 26)
ylabel('Stable P (kg/ha)', 'fontsize', 26)
lgd = legend('Reachable in 50 years','42 years','10 years','3 years','2 years','1 year', 'Target', 'Location', 'Northwest', 'fontsize', 26);
lgd.Color = 'white';
legend boxon

% Reachable sets with singular controls
figure
hold on
plot(Reach_array(4),'FaceColor',[0.4 0.4 0.7],'EdgeColor',[1 1 1],'FaceAlpha',1);
plot(Reach_array(3),'FaceColor',[0.2 0.2 0.5],'EdgeColor',[1 1 1],'FaceAlpha',1);
plot(Reach_array(2),'FaceColor',[0.1 0.1 0.3],'EdgeColor',[1 1 1],'FaceAlpha',1);
plot(Reach_array(1),'FaceColor',[0 0.0 0.1],'EdgeColor',[1 1 1],'FaceAlpha',1); 
plot(SafeInv,'FaceColor',[0.5 0.5 0.5],'FaceAlpha',1);
plot(Singular_array_all,'FaceColor',[0.98 0.85 0.54],'EdgeColor',[1 1 1],'FaceAlpha',1);
% plot initial condition estimate
plot([x1_init],[x2_init],'pentagram','MarkerFaceColor',[0.75 0 0.15],'MarkerEdgeColor','white', 'MarkerSize',50,'linewidth',2);
% plot formatting
xlim([mask_x 200])
ylim([350 1000])
ax = gca; % axes handle
ax.FontSize = 18; 
title('', 'fontsize', 26)
xlabel('Labile P (kg/ha)', 'fontsize', 26)
ylabel('Stable P (kg/ha)', 'fontsize', 26)
lgd = legend('Reachable in 4 years','3 years','2 years','1 year', 'Target', 'Singular regions','Location', 'Southwest', 'fontsize', 26);
lgd.Color = 'white';
legend boxoff

% Convex U calculation and plot lines
lgd.AutoUpdate = 'off';
c1 = 1 - alpha_3 - (beta_1*alpha_2 + beta_2*(1-alpha_2))*Q_bar;
c2 = alpha_3*alpha_1;
c3 = alpha_3;
c4 = 1 - alpha_3*alpha_1 - beta_2*Q_bar;

x1_min = (alpha_4*alpha_5 / (c1 - c2*c3/c4))^0.5 - alpha_5;

x1_1 = 130.4;  %148.5;
x1_2 = 104;  %126.1;
x1_3 = 61.29;  %90.42;
xint_1 = 930.2 + c3/c4 * x1_1; % 979.4
xint_2 = 945.6 + c3/c4 * x1_2; % 990.8
xint_3 = 973 + c3/c4 * x1_3; % 1012
hline = refline(-c3/c4,xint_1);
hline.Color = 'r';
hline.LineWidth = 2;
hline =refline(-c3/c4,xint_2);
hline.Color = 'r';
hline.LineWidth = 2;
hline =refline(-c3/c4,xint_3);
hline.Color = 'r';
hline.LineWidth = 2;

% Control Regions
% figure
% hold on
% plot(pmask,'FaceColor',[0 0.2 0.35],'FaceAlpha',1);
% plot(umaxback_plot,'FaceColor',[0.42 0.58 0.69],'FaceAlpha',1);
% plot(Singular_array_all,'FaceColor',[0.98 0.85 0.54],'FaceAlpha',1);
% plot(SafeInv,'FaceColor',[0.5 0.5 0.5],'FaceAlpha',1);
% % plot initial condition estimate
% plot([x1_init],[x2_init],'pentagram','MarkerFaceColor',[0.75 0 0.15],'MarkerEdgeColor','white', 'MarkerSize',25,'linewidth',2);
% % plot formatting
% xlim([30 220])
% ylim([mask_y 1800])
% ax = gca; % axes handle
% ax.FontSize = 18; 
% title('', 'fontsize', 26)
% xlabel('Labile P (kg/ha)', 'fontsize', 26)
% ylabel('Stable P (kg/ha)', 'fontsize', 26)
% lgd = legend('Minimum control region','Maximum control region','Singular region','Target', 'Location', 'Northwest', 'fontsize', 26);
% lgd.Color = 'white';
% legend boxon


%% Convex U calculation and plot lines
% c1 = 1 - alpha_3 - (beta_1*alpha_2 + beta_2*(1-alpha_2))*Q_bar;
% c2 = alpha_3*alpha_1;
% c3 = alpha_3;
% c4 = 1 - alpha_3*alpha_1 - beta_2*Q_bar;
% 
% x1_min = (alpha_4*alpha_5 / (c1 - c2*c3/c4))^0.5 - alpha_5;
% 
% x1_1 = 146.1;
% x1_2 = 111.9;
% x1_3 = 45.6;
% xint_1 = 1043 + c3/c4 * x1_1;
% xint_2 = 1065 + c3/c4 * x1_2;
% xint_3 = 1110 + c3/c4 * x1_3;
% hline = refline(-c3/c4,xint_1);
% hline.Color = 'r';
% hline =refline(-c3/c4,xint_2);
% hline.Color = 'r';
% hline =refline(-c3/c4,xint_3);
% hline.Color = 'r';
