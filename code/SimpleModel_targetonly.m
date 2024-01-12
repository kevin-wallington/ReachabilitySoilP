%% Mask for relevant region
mask_x = 30;
mask_y = 0;
mask_points = [mask_x mask_y;mask_x 15000;5000 15000;5000 mask_y];
pmask = polyshape(mask_points);

%% Create target set by back-propagating
% Safe and control invariant, i.e., can keep state in safe set
u_num = 20; %number of samples in control space
u_vec = linspace(u_min,u_max,u_num); 
horizon = 50; % number of seasons or years to simulate
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
SafeInv = Inv_array(50);

%% Plots
% % Safe, invariant, and recurrent sets
% figure
% hold on
% plot(SafeSet,'FaceColor',[0 0 0], 'FaceAlpha',1);
% plot(SafeInv,'FaceColor',[0.5 0.5 0.5], 'FaceAlpha',1);
% 
% xlim([0 450])
% ylim([000 1500])
% ax = gca; % axes handle
% ax.FontSize = 18; 
% title('', 'fontsize', 26)
% xlabel('Short-term P (kg/ha)', 'fontsize', 26)
% ylabel('Long-term P (kg/ha)', 'fontsize', 26)
% legend('Safe set','Target set - safe and invariant','Safe-Invariant-Recurrent set', 'Location', 'Northeast')
% legend boxoff