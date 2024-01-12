%% Mask - should be in place from prior routines...
% mask_x = 30;
% mask_points = [mask_x 0;mask_x 2500;400 2500;400 0];
% pmask = polyshape(mask_points);

%% Back propagate for u_max only
horizon = 30; % number of seasons or years to simulate
umax_region_array = repmat(polyshape, 1, horizon);
umin_region_array = repmat(polyshape, 1, horizon);
% pgon = SafeInv;
umax_pgon = SafeInv;
% umax_pgon_all = SafeInv;
umax_pgon_all = polyshape(); % debug, replaces above
umin_pgon = SafeInv;
% umin_pgon_all = SafeInv;
umin_pgon_all = polyshape(); % debug, replaces above
new_pgon = polyshape();
for jj=1:horizon % back propagation loop
    for ll=1:2
        samples = 0;
        if(ll==2)
            u=u_max;
            [cx,cy] = boundary(umax_pgon); % will wrap back to original pt to close polygon
        end
        if(ll==1)
            u=u_min;
            [cx,cy] = boundary(umin_pgon);
        end
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
        new_pgon=polyshape(u_poly_points);
        if(ll==2)
            if(jj>1); umax_pgon = subtract(new_pgon,umin_pgon_all); end % debug added if
            % new_pgon = subtract(new_pgon,umin_pgon);
            umax_pgon = subtract(umax_pgon,Singular_array_all);
            umax_pgon = intersect(umax_pgon,pmask);          
            umax_pgon_all = union(umax_pgon_all,umax_pgon);
            % umax_pgon = union(new_pgon,umax_pgon);
            % umax_pgon = subtract(umax_pgon,Singular_array_all);
            % umax_pgon = intersect(umax_pgon,pmask);
            umax_region_array(jj)=umax_pgon_all;
        end
        if(ll==1)
            umin_pgon = new_pgon;
            umin_pgon = subtract(umin_pgon,Singular_array_all);
            umin_pgon = intersect(umin_pgon,pmask);
            umin_pgon_all = union(umin_pgon_all,umin_pgon);
            % umin_pgon = union(new_pgon,umin_pgon);
            % umin_pgon = subtract(umin_pgon,Singular_array_all);
            % umin_pgon = intersect(umin_pgon,pmask);
        end
    end
end

umax_region_pgon = umax_region_array(30);

%% 
background_points = [mask_x mask_y;mask_x 600;100 600;100 mask_y];
pgon = polyshape(background_points);
umax_region_pgon = union(umax_region_pgon, pgon); % debug, new
%% Plot
figure
hold on
plot(pmask,'FaceColor',[0 0.2 0.35],'FaceAlpha',1);
plot(umax_region_pgon,'FaceColor',[0.42 0.58 0.69],'FaceAlpha',1);
plot(Singular_array_all,'FaceColor',[0.98 0.85 0.54],'FaceAlpha',1);
plot(SafeInv,'FaceColor',[0.5 0.5 0.5],'FaceAlpha',1);
% plot initial condition estimate
plot([x1_init],[x2_init],'pentagram','MarkerFaceColor',[0.75 0 0.15],'MarkerEdgeColor','white', 'MarkerSize',50,'linewidth',2);
% plot formatting
xlim([30 220])
ylim([mask_y 1800])
ax = gca; % axes handle
ax.FontSize = 18; 
title('', 'fontsize', 26)
xlabel('Labile P (kg/ha)', 'fontsize', 26)
ylabel('Stable P (kg/ha)', 'fontsize', 26)
lgd = legend('Minimum control region','Maximum control region','Singular region','Target', 'Location', 'Northwest', 'fontsize', 26);
lgd.Color = 'white';
legend boxon