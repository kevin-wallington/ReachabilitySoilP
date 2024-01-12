%% Mask - should be in place from prior routines...
% mask_x = 30;
mask_points = [mask_x 0;mask_x 15000;5000 15000;5000 0];
pmask = polyshape(mask_points);
% mask_points_maxback = [mask_x mask_y;mask_x 15000;x1_min+50 15000;x1_min+50 mask_y];% DEBUG 07/20/23 
pmask = polyshape(mask_points);
% maxback_mask = polyshape(mask_points_maxback); % DEBUG 07/20/23

%% Back propagate for u_max only
horizon = 20; % number of seasons or years to simulate - reduced to 20 for sensitivity purposes, 07/21/23
umax_region_array = repmat(polyshape, 1, horizon);
umin_region_array = repmat(polyshape, 1, horizon);
% pgon = SafeInv;
umax_pgon = SafeInv;
% [cx,cy] = boundary(umax_pgon);
% max_x = max(cx);
% mask_points = [mask_x 0;mask_x 15000;max_x 15000;max_x 0];
% pmask = polyshape(mask_points);
% umax_pgon_all = SafeInv; %%% DEBUG 07/20/22
umax_pgon_all = polyshape();
% new_pgon = polyshape();
u = u_max;
for jj=1:horizon % back propagation loop
%     [cx,cy] = boundary(umax_pgon_all); %%% DEBUG 07/20/23
    [cx,cy] = boundary(umax_pgon); % will wrap back to original pt to close polygon
    poly_points=[];
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
    u_poly_points = zeros(samples,2);
    eliminate = [];
    for ii=1:samples
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
            % max_minx1 = max(max_minx1, min(root_out));
            new_x2 = poly_points(ii,2) / A(2,2) -  (A(2,1)/A(2,2))*new_x1;
            u_poly_points(ii,:) = [new_x1 new_x2]; 
        end
    end
    u_poly_points(eliminate,:) = [];
    u_poly_points(abs(u_poly_points)<0.1)=0;
    umax_pgon=polyshape(u_poly_points);
%     umax_pgon = subtract(new_pgon,umin_pgon_all);
%     umax_pgon = subtract(umax_pgon,Singular_array_all);
    umax_pgon = intersect(umax_pgon,pmask);      
    umax_pgon_all = union(umax_pgon_all,umax_pgon);
%     umax_pgon_all = intersect(umax_pgon_all,maxback_mask);   % DEBUG 07/20/23
    umax_region_array(jj)=umax_pgon_all;
%     if jj<=4 %%%%COMMENT OUT??? DEBUG 07/20/23
%         umax_pgon_all = union(umax_pgon_all, Singular_array(jj));
%     end
%     umax_pgon_all = intersect(umax_pgon_all,maxback_mask);   % DEBUG 07/20/23
%     umax_pgon_all = union(umax_pgon_all,SafeInv);   % DEBUG 07/20/23   
end

umax_region_pgon = umax_region_array(20); % reduced to 20 for sensitivity purposes, 07/21/23

%% Plot
% figure
% hold on
% plot(umax_region_pgon,'FaceColor',[0 0.2 0.35],'FaceAlpha',1);
% plot(SafeInv,'FaceColor',[0.5 0.5 0.5], 'FaceAlpha',1); 
% plot(Singular_array_all,'FaceColor',[0.98 0.85 0.54],'FaceAlpha',1);
% % plot formatting
% xlim([30 250])
% ylim([0 1500])
% ax = gca; % axes handle
% ax.FontSize = 18; 
% title('', 'fontsize', 26)
% xlabel('Short-term P (kg/ha)', 'fontsize', 26)
% ylabel('Long-term P (kg/ha)', 'fontsize', 26)