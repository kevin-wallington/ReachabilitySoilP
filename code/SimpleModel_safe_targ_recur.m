%% Create safe set
% 97.5% chance of satisfying water quality, satisfy nutrient deficit constraint
Def_lim = 0.75; % nutrient deficit limit (i.e. 75% of max)
Load_lim = 0.55*Q_bar*((beta_1*alpha_2+beta_2*(1-alpha_2))*x1_init + beta_2*x2_init); % 226 and 1850 are estimated IC had 0.6 *
% Load_lim = 1.08; % maximum P losses, 40% reduction from 1.8 kg/ha ; previously had 1.75 kg/ha
% x3_nom = 0.6; %0.5; 
% ps_nom = 0.0; %0.4;
x1_min = Def_lim/(1-Def_lim)*alpha_5;
x2_max = (Load_lim - Q_bar*(beta_1*alpha_2+beta_2*(1-alpha_2))*x1_min) / (Q_bar*beta_2);
x1_max = Load_lim /(Q_bar*(beta_1*alpha_2+beta_2*(1-alpha_2)));
cx = [x1_min, x1_min, x1_max]; % corner points, x-coord
cy = [0, x2_max, 0]; % corner points, y-coord
poly_points = [cx; cy]';
pgon=polyshape(poly_points);
SafeSet=pgon;

%% Masks
mask_x = 30;
mask_y = 0;
mask_points = [mask_x mask_y;mask_x 2500;1500 2500;1500 mask_y];
pmask = polyshape(mask_points);

%% Create target set by back-propagating
% Safe and control invariant, i.e., can keep state in safe set 
horizon = 100; % number of seasons or years to simulate
Inv_array = repmat(polyshape, 1, horizon);
pgon = SafeSet;
for jj=1:horizon % back propagation loop
    poly_points=[];
    [cx,cy] = boundary(pgon); % will wrap back to original pt to close polygon
    num_pts = length(cx);
    for kk=1:(num_pts-1) %generate points along boundary of polygon
       d=pdist([cx(kk) cy(kk); cx(kk+1) cy(kk+1)]);
       seg_pts = floor(d/5)+2;
       x_seg = linspace(cx(kk),cx(kk+1),seg_pts);
       y_seg = linspace(cy(kk),cy(kk+1),seg_pts);
       x_seg(end)=[]; y_seg(end)=[];
       new_pts = [x_seg; y_seg];
       poly_points = [poly_points' new_pts]';
    end    
    samples=length(poly_points);
    poly_points(abs(poly_points)<0.1)=0; % otherwise, points very close to zero
    umax_poly_points = zeros(samples,2);
    umin_poly_points = zeros(samples,2);
    for ll=1:2
        eliminate = [];
        if(ll==1); u=u_max; end
        if(ll==2); u=u_min; end        
        for ii=1:samples % back propagate each point
            c5 = A(1,1) - A(1,2) * A(2,1) / A(2,2);
            c6 = (A(1,2) / A(2,2)) * poly_points(ii,2) + u;
            a = c5; 
            b = c5*alpha_5 - alpha_4 + c6 - poly_points(ii,1);
            c = c6*alpha_5 - alpha_5*poly_points(ii,1);
            root_out = roots([a b c]);
            if (~isreal(root_out))
                eliminate = [eliminate, ii]; 
            else
                new_x1 = max(root_out);
            end
            new_x2 = poly_points(ii,2) / A(2,2) -  (A(2,1)/A(2,2))*new_x1;
            if(ll==1); umax_poly_points(ii,:) = [new_x1 new_x2]; end
            if(ll==2); umin_poly_points(ii,:) = [new_x1 new_x2]; end 
        end
        if(ll==1); umax_poly_points(eliminate,:) = []; end
        if(ll==2); umin_poly_points(eliminate,:) = []; end
    end
    umax_poly_points(abs(umax_poly_points)<0.1)=0;
    umin_poly_points(abs(umin_poly_points)<0.1)=0;
    umax_pgon=polyshape(umax_poly_points);
    umin_pgon=polyshape(umin_poly_points);
    new_pgon = union(umax_pgon,umin_pgon);
    new_pgon = intersect(new_pgon,pmask);
    [Rch_cx, Rch_cy] = convU(new_pgon, A, mask_x, mask_y);  
    pgon = polyshape(Rch_cx, Rch_cy);
    pgon = intersect(pgon,SafeSet); % only want priors that are safe (and map to safe)
    Inv_array(jj)=pgon;
end
SafeInv = Inv_array(100);

%% Create recurrent set by forward propogating
u_num = 20; %number of samples in control space
u_vec = linspace(u_min,u_max,u_num);

% Forward propagate
horizon = 200; % number of seasons or years to simulate - CONVERGES VERY SLOWLY
Recur_array = repmat(polyshape, 1, horizon);
pgon = SafeInv;
for jj=1:horizon % back propagation loop
    poly_points=[];
    [cx,cy] = boundary(pgon); % will wrap back to original pt to close polygon
    num_pts = length(cx);
    for kk=1:(num_pts-1) %generate points along boundary of polygon
       d=pdist([cx(kk) cy(kk); cx(kk+1) cy(kk+1)]);
       seg_pts = floor(d/5)+2;
       x_seg = linspace(cx(kk),cx(kk+1),seg_pts);
       y_seg = linspace(cy(kk),cy(kk+1),seg_pts);
       x_seg(end)=[]; y_seg(end)=[];
       new_pts = [x_seg; y_seg];
       poly_points = [poly_points' new_pts]';
    end    
    samples=length(poly_points);
    poly_points(abs(poly_points)<0.1)=0; % otherwise, points very close to zero
    umax_poly_points = zeros(samples,2);
    umin_poly_points = zeros(samples,2);
    for ll=1:2%:u_num
        if(ll==1); u=u_max; end
        if(ll==2); u=u_min; end 
        for ii=1:samples 
            f_nonlin = [(-alpha_4*(poly_points(ii,1)/(poly_points(ii,1)+alpha_5))); 0];
            b = [u;0]; 
            if(ll==1)
                umax_poly_points(ii,:) = (A * poly_points(ii,:)' + f_nonlin + b)'; % update state
            end
            if(ll==2)
                umin_poly_points(ii,:) = (A * poly_points(ii,:)' + f_nonlin + b)'; % update state
            end
        end
    end
    umax_poly_points(abs(umax_poly_points)<0.1)=0;
    umin_poly_points(abs(umin_poly_points)<0.1)=0;
    umax_pgon=polyshape(umax_poly_points);
    umin_pgon=polyshape(umin_poly_points);
    new_pgon = union(umax_pgon,umin_pgon);
    [Rch_cx, Rch_cy] = convU_forward(new_pgon, A, mask_x, mask_y); 
    pgon = polyshape(Rch_cx, Rch_cy);
    pgon = intersect(pgon,SafeInv); % only want priors that are safe (and map to safe)
    Recur_array(jj)=pgon;
end
SafeRecur = Recur_array(200);

%% Identify control values that render invariant set as such
% Finer control resolution
u_num = 500; %number of samples in control space
u_vec = linspace(u_min,u_max,u_num); 
% Sample boundary of SafeInv set
pgon=SafeInv; % Safe_Inv created by SimpleModel_Invariance.m
poly_points=[];
[cx,cy] = boundary(pgon); % will wrap back to original pt to close polygon
num_pts = length(cx);
for kk=1:(num_pts-1) %generate points along boundary of polygon
   d=pdist([cx(kk) cy(kk); cx(kk+1) cy(kk+1)]);
   seg_pts = floor(d)+2; % or floor(d/5)
   x_seg = linspace(cx(kk),cx(kk+1),seg_pts);
   y_seg = linspace(cy(kk),cy(kk+1),seg_pts);
   x_seg(end)=[]; y_seg(end)=[];
   new_pts = [x_seg; y_seg];
   poly_points = [poly_points' new_pts]';
end    
samples=length(poly_points);
poly_points(abs(poly_points)<0.1)=0; % otherwise, points very close to zero

% Forward propagate one step
pgon=SafeInv;
polyout = polybuffer(pgon,0.04); % create tolerance for numerical error
[cx,cy] = boundary(polyout); 
horizon = 1; % number of seasons or years to simulate - CONVERGES VERY SLOWLY
InvCont_array = repmat(polyshape, 1, u_num);
u_poly_points = cell(u_num,1);
u_poly_points(:) = {zeros(samples, 2)}; 
u_in_on = cell(u_num,1);
u_in_on(:) = {zeros(samples, 2)}; 
u_bounds = zeros(samples, 2);  % first column is minimum u which satisfies, second is max
for ll=1:u_num
    u = u_vec(ll);
    for ii=1:samples 
        f_nonlin = [(-alpha_4*(poly_points(ii,1)/(poly_points(ii,1)+alpha_5))); 0];
        b = [u;0]; 
        u_poly_points{ll}(ii,:) = (A * poly_points(ii,:)' + f_nonlin + b)'; % update state
    end
    u_poly_points{ll}(abs(u_poly_points{ll})<0.1)=0;
    [u_in_on{ll}(:,1),u_in_on{ll}(:,2)] = inpolygon(u_poly_points{ll}(:,1),u_poly_points{ll}(:,2),cx,cy); % first column returns 1 if point is in poly, second returns 1 if on boundary
    for ii=1:samples 
        if (u_in_on{ll}(ii,1)==1 || u_in_on{ll}(ii,2)==1) 
            u_bounds(ii,2)=u;
            if (u_bounds(ii,1)==0)
               u_bounds(ii,1)=u;
            end
        end
    end
end

%% Plots
% % Safe, invariant, and recurrent sets
figure
hold on
plot(SafeSet,'FaceColor',[0 0 0], 'FaceAlpha',1);
plot(SafeInv,'FaceColor',[0.5 0.5 0.5], 'FaceAlpha',1);
plot(SafeRecur,'FaceColor',[0.98 0.85 0.54], 'FaceAlpha',1);

xlim([100 265])
ylim([000 1000])
ax = gca; % axes handle
ax.FontSize = 18; 
title('', 'fontsize', 26)
xlabel('Labile P (kg/ha)', 'fontsize', 26)
ylabel('Stable P (kg/ha)', 'fontsize', 26)
legend('Safe set','Target set - safe and invariant','Safe-Invariant-Recurrent set', 'Location', 'Northeast', 'fontsize', 26)
legend boxoff
% 
% % plot equilibrium points (for some control level)
% x1_eq_min = 120; % 120 for base
% x1_eq_max = 130; % 130 for base
% num = 11;
% x1_range = linspace(x1_eq_min, x1_eq_max,num);
% x2_eq_range = alpha_3*x1_range/(beta_2*Q_bar+alpha_3*alpha_1);
% pt_min = [x1_range(1) , x2_eq_range(1)];
% pt_max = [x1_range(num) , x2_eq_range(num)];
% f_nonlin = [(-alpha_4*pt_min(1))/(pt_min(1)+alpha_5); 0];
% ptprop = (A * pt_min' + f_nonlin)';
% eq_umin = x1_eq_min - ptprop(1) - 7;
% f_nonlin = [(-alpha_4*pt_max(1))/(pt_max(1)+alpha_5); 0];
% ptprop = (A * pt_max' + f_nonlin)';
% eq_umax = x1_eq_max - ptprop(1) - 7;
% plot(x1_range,x2_eq_range,'k');
% 
% % Controls along invariant boundary
figure
hold on
scatter(poly_points(:,1), poly_points(:,2), [], u_bounds(:,1), 'Filled'); % U_min
caxis([u_min, u_max]);
xlim([100 220])
ylim([300 1000])
ax = gca; % axes handle
ax.FontSize = 24; 
title('', 'fontsize', 32)
xlabel('Labile P (kg/ha)', 'fontsize', 32)
ylabel('Stable P (kg/ha)', 'fontsize', 32)
a = colorbar;
a.Label.String = 'Minimum fertilizer rate for invariance';
a.FontSize = 24;

figure
hold on
scatter(poly_points(:,1), poly_points(:,2), [], u_bounds(:,2), 'Filled'); % U_max
caxis([u_min, u_max]);
xlim([100 220])
ylim([300 1000])
ax = gca; % axes handle
ax.FontSize = 24; 
title('', 'fontsize', 26)
xlabel('Labile P (kg/ha)', 'fontsize', 32)
ylabel('Stable P (kg/ha)', 'fontsize', 32)
a = colorbar;
a.Label.String = 'Maximum fertilizer rate for invariance';
a.FontSize = 24;