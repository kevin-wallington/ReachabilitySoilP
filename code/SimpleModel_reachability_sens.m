%% Create target polygon points vector and mask for state space
pgon= SafeInv;
umax_pgon = SafeInv;
umaxback_pgon = SafeInv;
[cx,cy] = boundary(pgon);
mask_y= min(cy); % DEBUG 07/20/23 - cuts off bottom of polygons so don't have issues
if(mask_y<=0.1);mask_y = 0; end
if(mask_x<=0.1);mask_x = 0; end
% mask_x = 30;
% mask_y = 0;
mask_points = [mask_x mask_y;mask_x 15000;5000 15000;5000 mask_y];
% mask_points_maxback = [mask_x mask_y;mask_x x2_max;x1_min x2_max;x1_min mask_y]; % DEBUG for l=3
pmask = polyshape(mask_points);
% maxback_mask = polyshape(mask_points_maxback); % DEBUG for l=3

%% Back propagate
horizon = 20; % number of seasons or years to simulate
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
    premax_poly_points=[]; % DEBUG umax_pgon instead of umaxback_pgon , 07/20/23
    [mbx,mby] = boundary(umax_pgon); % DEBUG umax_pgon instead of umaxback_pgon , 07/20/23
%     maxback_poly_points=[];
%     [mbx,mby] = boundary(umaxback_pgon); % will wrap back to original pt to close polygon
    num_pts = length(mbx);
    for kk=1:(num_pts-1) %generate points along boundary of polygon
       d=pdist([mbx(kk) mby(kk); mbx(kk+1) mby(kk+1)]);
       seg_pts = floor(d/10)+2;
       x_seg = linspace(mbx(kk),mbx(kk+1),seg_pts);
       y_seg = linspace(mby(kk),mby(kk+1),seg_pts);
       x_seg(end)=[]; y_seg(end)=[];
       new_pts = [x_seg; y_seg];
       premax_poly_points = [premax_poly_points' new_pts]';  % DEBUG umax_pgon instead of umaxback_pgon , 07/20/2
%        maxback_poly_points = [maxback_poly_points' new_pts]';
    end    
    premax_poly_points(abs(premax_poly_points)<0.1)=0; % DEBUG umax_pgon instead of umaxback_pgon , 07/20/2
%     maxback_poly_points(abs(maxback_poly_points)<0.1)=0; % otherwise, points very close to zero
    % Backpropogate with bang-bang control
    umax_poly_points = zeros(samples,2);
    umin_poly_points = zeros(samples,2);
%     umax3_poly_points = zeros(samples,2); % for l=3
    %%% DEBUG adding l=3
    for ll=1:2 % 3
        %%%%% DEBUG - original has linearization not elimination
        eliminate = [];
        poly_points=[];
        if(ll==1)
            u=u_max; 
            poly_points=premax_poly_points;  % DEBUG umax_pgon instead of umaxback_pgon , 07/20/2
%             poly_points=maxback_poly_points;  % Replace with full set for
            % sensitivity analysis so code runs without errors when two
            % polygons
%             poly_points=other_poly_points;
        end
        if(ll==2)
            u=u_min; 
            poly_points=other_poly_points; 
        end
%         if(ll==3) %%%%DEBUG
%             u=u_max; 
%             poly_points=maxback_poly_points; 
%         end
        samples=length(poly_points);
        for ii=1:samples % back propagate each point
            c5 = A(1,1) - A(1,2) * A(2,1) / A(2,2);
            c6 = (A(1,2) / A(2,2)) * poly_points(ii,2) + u;
            a = c5; 
            b = c5*alpha_5 - alpha_4 + c6 - poly_points(ii,1);
            c = c6*alpha_5 - alpha_5*poly_points(ii,1);
            root_out = roots([a b c]);
            if (~isreal(root_out))
                %%%%% DEBUG - original has linearization, not elimination
                eliminate = [eliminate, ii]; 
%                 pup = alpha_4*poly_points(ii,1)/(poly_points(ii,1)+alpha_5);
%                 d_pup = alpha_4*alpha_5/(poly_points(ii,1)+alpha_5)^2;
%                 if(d_pup>=0.95*A(1,1));d_pup=0.95*A(1,1);end % to prevent instability near x1 = 0
%                 b1 = u - pup + d_pup*poly_points(ii,1);
%                 b = [b1;0]; 
%                 back_prop = poly_points(ii,:)' - b;
%                 A_loc = A;
%                 A_loc(1,1)=A_loc(1,1)-d_pup;
%                 init_est = (A_loc\back_prop)';
%                 new_x1 = init_est(1);
            else
                new_x1 = max(root_out);
                % max_minx1 = max(max_minx1, min(root_out));
            end
            new_x2 = poly_points(ii,2) / A(2,2) -  (A(2,1)/A(2,2))*new_x1;
            if(ll==1); umax_poly_points(ii,:) = [new_x1 new_x2]; end
            if(ll==2); umin_poly_points(ii,:) = [new_x1 new_x2]; end 
%             if(ll==3); umax3_poly_points(ii,:) = [new_x1 new_x2]; end 
        end
        % DEBUG - original has linearization not elimination
        if(ll==1); umax_poly_points(eliminate,:) = []; end
        if(ll==2); umin_poly_points(eliminate,:) = []; end
%         if(ll==3); umax3_poly_points(eliminate,:) = []; end
    end
    umax_poly_points(abs(umax_poly_points)<0.1)=0;
    umin_poly_points(abs(umin_poly_points)<0.1)=0;
%     umax3_poly_points(abs(umax3_poly_points)<0.1)=0; % for l=3
    umax_pgon=polyshape(umax_poly_points);
    umin_pgon=polyshape(umin_poly_points);
%     umax3_pgon=polyshape(umax3_poly_points); % for l=3
    new_pgon = union(umax_pgon, umin_pgon); % NEW
%     new_pgon = union(pgon,umax_pgon);
%     new_pgon = union(new_pgon,umin_pgon);
    new_pgon = intersect(new_pgon,pmask);
    umax_pgon = intersect(umax_pgon,pmask); % DEBUG 07/20/23, just use umax_pgon for backprop
%     umaxback_pgon = union(umax_pgon,SafeInv); % DEBUG - changed to max3. DEBUG AGAIN 07/20/23, just use umax_pgon for backprop
%     umaxback_pgon = intersect(umaxback_pgon,pmask); %
    BBreach_array(jj)=new_pgon;
    [Rch_cx, Rch_cy] = convU(BBreach_array(jj), A, mask_x, mask_y);  
%     pgon = polyshape(Rch_cx, Rch_cy); 
    new_pgon = polyshape(Rch_cx, Rch_cy); % NEW
    pgon = union(pgon,new_pgon); %NEW
    Reach_array(jj)=pgon;
    Singular_array(jj) = subtract(Reach_array(jj),BBreach_array(jj));
    Singular_array_control(jj) = Singular_array(jj);
    Singular_array(jj) = subtract(Singular_array(jj),SafeInv);
    Singular_array_control(jj) = Singular_array(jj);
    if (jj>1);Singular_array(jj) = subtract(Singular_array(jj),Reach_array(jj-1));end % new on 07/21/23
%     extra_sing = subtract(Reach_array(jj),umaxback_pgon); % DEBUG for l=3
%     extra_sing = intersect(extra_sing,maxback_mask); % DEBUG for l=3
%     Singular_array(jj) = union(Singular_array(jj),extra_sing); % DEBUG for l=3
end

%% Singular region excludes bang-bang reachable region
Singular_array_all = Singular_array(1);
for ii=2:8 % small regions due to numerical errors at 7 steps (4 is enough for meaningful regions in standard params) 
    % Check results with 0 hit time (in optimal scenario) with more or less
    % singular array steps (e.g. 2:8)
    Singular_array_all = union(Singular_array_all,Singular_array(ii));
end
Singular_array_all = subtract(Singular_array_all,SafeInv);