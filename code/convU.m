function [Rch_cx, Rch_cy] = convU(BBpgon, A, mask_x, mask_y)
% Remove points between "tips" of u_max and u_min regions from polyshape because intermediated controls can reach from that space (this is taking convex hull of just that one region)
% Note: this should really apply at bottom tips of reach sets too, butnot worth coding because those aren't realistic part of state space and won't even appear in figures
minx_maxy = 5000;
maxx_maxy = 0;
[BBcx,BBcy] = boundary(BBpgon);
[BBcx,BBcy] = poly2ccw(BBcx,BBcy); %counterclockwise
BBcx(abs(BBcx)<0.1)=0;
BBcy(abs(BBcy)<0.1)=0;
% maxy = max(BBcy); % identify maximum y coordinate - the tip
% miny = min(BBcy); % identify maximum y coordinate - the tip
% k_max = find(BBcy==maxy);
% Find polygon corner point closest to level line (for u convexity)
p0 = [BBcx, BBcy];
p1 = [-200, (0-(A(2,1)/A(2,2))*(0-200))]; % arbitrary point on line for creating vector
p2 = [800, (0-(A(2,1)/A(2,2))*(0+800))]; % arbitrary point on line for creating vector
d = abs(((p2(:,1)-p1(:,1)).*(p1(:,2)-p0(:,2))-(p2(:,2)-p1(:,2)).*(p1(:,1)-p0(:,1))) ./ sqrt((p2(:,1)-p1(:,1)).^2+(p2(:,2)-p1(:,2)).^2));
d_max = max(d); % identify maximum distance
d_min = min(d); 
k_max = find(ismembertol(d(2:end),d_max,0.000001)==1)'+1; % index of pt with max distance, numerical tolerance
k_min = find(ismembertol(d(2:end),d_min,0.000001)==1)'+1; 
if (sum(isnan(BBcx))==0) % only one polygon present
    for ii=k_max
        if (BBcx(ii)<minx_maxy)
            minx_maxy=BBcx(ii);
            k_topleft = ii;
        end
    end
    if k_topleft < length(BBcx) % that is, does not begin and end with maximum y point
%     if length(k_max)==1 % that is, does not begin and end with maximum y point
        BBcx = [BBcx(k_topleft:end); BBcx(2:k_topleft)]; % change so that do begin and end with maximum y point, used to have 'k_max' not 'k_topleft'
        BBcy = [BBcy(k_topleft:end); BBcy(2:k_topleft)];
    end
    % Find polygon corner point closest to level line (for u convexity)
    p0 = [BBcx, BBcy];
    p1 = [-200, (0-(A(2,1)/A(2,2))*(0-200))]; % arbitrary point on line for creating vector
    p2 = [800, (0-(A(2,1)/A(2,2))*(0+800))]; % arbitrary point on line for creating vector
    d = abs(((p2(:,1)-p1(:,1)).*(p1(:,2)-p0(:,2))-(p2(:,2)-p1(:,2)).*(p1(:,1)-p0(:,1))) ./ sqrt((p2(:,1)-p1(:,1)).^2+(p2(:,2)-p1(:,2)).^2));
%     d_max = max(d); % identify maximum distance
    k_max = find(ismembertol(d(2:end),d_max,0.000001)==1)'+1; % index of pt with max distance, numerical tolerance
    k_min = find(ismembertol(d(2:end),d_min,0.000001)==1)'+1; % index of pt with max distance, numerical tolerance
    if length(k_max)>=2 % NEW 07/25/23 changed from ==2 to >= 2, for degenerate case where extra point within tolerance
        % eliminate point
        k_max_2 = k_max(end); % max location should now be at beginning and end, NEW on 07/25/23, changed k_max(2) to k_max(:), for degenerate case where extra point within tolerance
        k_max_1 = k_max(1); % the k_max_end point is not included so largest point is prior from ccw direction
        Rch_cx = BBcx;
        Rch_cy = BBcy;
        Rch_cx(k_max_1+1:k_max_2-1)=[]; % remove points that are not part of convex hull, just for this specific region
        Rch_cy(k_max_1+1:k_max_2-1)=[]; 
    end
    if length(k_max)==1 && ~ismembertol(BBcx(k_max),mask_x,0.000001)
        % add x_mask intercept as polygon point
        k_xmask = find(ismembertol(BBcx,mask_x,0.000001)==1)';
        y_max = 0;
        for ii=k_xmask
            if BBcy(ii)>y_max
               y_max = BBcy(ii);
            end
        end
        BBcx_alt = BBcx; BBcy_alt = BBcy;
        k_maxmask = k_xmask(find(BBcy(k_xmask)==y_max))';
        for ii=k_maxmask
           BBcy_alt(ii)=(BBcy(k_max)-(A(2,1)/A(2,2))*(mask_x-BBcx(k_max))); % where c3/c4 line intersects
        end
        if length(k_maxmask)==1 % i.e., therefore the max y point and first entry for ccw rotation is from u_min
            k_max_2 = k_maxmask;
            k_max_1 = 1; % max location should now be at beginning and end
        else
            k_max_2 = length(BBcx_alt); % max location should now be at beginning and end
            k_max_1 = k_max(1);
        end
        Rch_cx = BBcx_alt;
        Rch_cy = BBcy_alt;
        Rch_cx(k_max_1+1:k_max_2-1)=[]; % remove points that are not part of convex hull, just for this specific region
        Rch_cy(k_max_1+1:k_max_2-1)=[]; 
    end
    if length(k_max)==1 && ismembertol(BBcx(k_max),mask_x,0.000001)
        %remove all points b/w max x2 on mask and min x2 on mask (moving CCW)
        k_xmask = find(ismembertol(BBcx,mask_x,0.000001)==1)';
        y_max = 0;
        y_min = 15000;
        for ii=k_xmask
            if BBcy(ii)>y_max
               y_max = BBcy(ii);
               k_maxmask = ii; 
            end
            if BBcy(ii)<y_min
               y_min = BBcy(ii);
               k_minmask = ii; 
            end
        end
        Rch_cx = BBcx;
        Rch_cy = BBcy;
        Rch_cx(k_maxmask+1:k_minmask-1)=[]; % remove points that are not part of convex hull, just for this specific region
        Rch_cy(k_maxmask+1:k_minmask-1)=[];
    end  
    % Same treatment for bottom of polygon
    p0 = [Rch_cx, Rch_cy];
    d = abs(((p2(:,1)-p1(:,1)).*(p1(:,2)-p0(:,2))-(p2(:,2)-p1(:,2)).*(p1(:,1)-p0(:,1))) ./ sqrt((p2(:,1)-p1(:,1)).^2+(p2(:,2)-p1(:,2)).^2));
    d_min = min(d); % identify minimum distance
    k_min = find(ismembertol(d(2:end),d_min,0.000001)==1)'+1; % index of pt with max distance, numerical tolerance
    if length(k_min)>=2  % NEW 07/25/23 changed from ==2 to >= 2, for degenerate case where extra point within tolerance
        % eliminate point
        k_min_2 = k_min(:); % max location should now be at beginning and end,NEW 07/25/23 changed from k_min(2) to k_min(end), for degenerate case where extra point within tolerance
        k_min_1 = k_min(1); % the k_max_end point is not included so largest point is prior from ccw direction
        Rch_cx(k_min_1+1:k_min_2-1)=[]; % remove points that are not part of convex hull, just for this specific region
        Rch_cy(k_min_1+1:k_min_2-1)=[]; 
    end
    if length(k_min)==1 && ~ismembertol(Rch_cy(k_min),mask_y,0.000001) 
        % add y_mask intercept as polygon point
        k_ymask = find(ismembertol(Rch_cy,mask_y,0.000001)==1)';
        x_min = 5000;
        for ii=k_ymask 
            if Rch_cx(ii)<x_min
               x_min = Rch_cx(ii);
               k_minmask = ii;
            end
        end
        Rchcx_alt = Rch_cx; Rchcy_alt = Rch_cy;
%         k_minmask = k_ymask(find(Rch_cx(k_ymask)==x_min))';
        for ii=k_minmask
           Rchcx_alt(ii)=(Rch_cx(k_min)-(A(2,2)/A(2,1))*(mask_y-Rch_cy(k_min))); % where c3/c4 line intersects
        end
        k_min_2 = k_minmask;
        k_min_1 = k_min(1); % max location should now be at beginning and end
        Rch_cx = Rchcx_alt;
        Rch_cy = Rchcy_alt;
        Rch_cx(k_min_1+1:k_min_2-1)=[]; % remove points that are not part of convex hull, just for this specific region
        Rch_cy(k_min_1+1:k_min_2-1)=[]; 
    end
    if length(k_min)==1 && ismembertol(Rch_cy(k_min),mask_y,0.000001)
        %remove all points b/w max x2 on mask and min x2 on mask (moving CCW)
        k_ymask = find(ismembertol(Rch_cy,mask_y,0.000001)==1)';
        if (k_ymask(end)-k_ymask(1))>1
            x_max = 0;
            x_min = 5000;
            for ii=k_ymask
                if Rch_cx(ii)<x_min
                   x_min = Rch_cx(ii);
                   k_minmask = ii; %??
                end
                if Rch_cx(ii)>x_max
                   x_max = Rch_cx(ii);
                   k_maxmask = ii; %??
                end
            end
%             Rch_cx = BBcx;
%             Rch_cy = BBcy;
            Rch_cx(k_minmask+1:k_maxmask-1)=[]; % remove points that are not part of convex hull, just for this specific region
            Rch_cy(k_minmask+1:k_maxmask-1)=[]; 
        end
    end 
else % two polygons
    nan_loc = find(isnan(BBcx));
    for ii=k_max
        if (BBcx(ii)<minx_maxy)
            minx_maxy=BBcx(ii);
            k_topleft = ii;
        end
        if (BBcx(ii)>maxx_maxy)
            maxx_maxy=BBcx(ii);
            k_topright = ii;
        end
    end
    if k_topleft==k_topright % just one point with max y - will be in right poly
        if k_topright>nan_loc
            maxy1 = max(d(1:(nan_loc-1))); % identify maximum distance - the tip
            k_max1 = find(d(2:(nan_loc-1))==maxy1)+1;
            k_topleft = k_max1(1); 
        else
            maxy1 = max(d((nan_loc+1):end)); 
            k_max1 = find(d((nan_loc+2):end)==maxy1)+nan_loc+1;
            k_topleft = k_max1(1); 
        end
    end  
    if k_topleft<nan_loc
        if ~(k_topleft==1) % that is, does not begin and end with maximum y point
            BBcx1 = [BBcx(k_topleft:(nan_loc-1)); BBcx(2:k_topleft)]; % change so that do begin and end with maximum y point
            BBcy1 = [BBcy(k_topleft:(nan_loc-1)); BBcy(2:k_topleft)];
        else
            BBcx1 = BBcx(1:(nan_loc-1));
            BBcy1 = BBcy(1:(nan_loc-1));
        end
        if ~(k_topright==(nan_loc+1)) % that is, does not begin and end with maximum y point
            BBcx2 = [BBcx(k_topright:end); BBcx((nan_loc+2):k_topright)]; % change so that do begin and end with maximum y point
            BBcy2 = [BBcy(k_topright:end); BBcy((nan_loc+2):k_topright)];
        else
            BBcx2 = BBcx((nan_loc+1):end);
            BBcy2 = BBcy((nan_loc+1):end);
        end
    else
        if ~(k_topleft==(nan_loc+1)) % that is, does not begin and end with maximum y point
            BBcx1 = [BBcx(k_topleft:end); BBcx((nan_loc+2):k_topleft)]; % change so that do begin and end with maximum y point
            BBcy1 = [BBcy(k_topleft:end); BBcy((nan_loc+2):k_topleft)];
        else
            BBcx1 = BBcx((nan_loc+1):end);
            BBcy1 = BBcy((nan_loc+1):end);
        end
        if ~(k_topright==1) % that is, does not begin and end with maximum y point
            BBcx2 = [BBcx(k_topright:(nan_loc-1)); BBcx(2:k_topright)]; % change so that do begin and end with maximum y point
            BBcy2 = [BBcy(k_topright:(nan_loc-1)); BBcy(2:k_topright)];
        else
            BBcx2 = BBcx(1:(nan_loc-1));
            BBcy2 = BBcy(1:(nan_loc-1));
        end  
    end
    BBcx = [BBcx1; NaN; BBcx2];
    BBcy = [BBcy1; NaN; BBcy2];
    nan_loc = find(isnan(BBcx));
    % Find polygon corner point closest to level line (for u convexity)
    p0 = [BBcx, BBcy];
    p1 = [-200, (0-(A(2,1)/A(2,2))*(0-200))]; % arbitrary point on line for creating vector
    p2 = [800, (0-(A(2,1)/A(2,2))*(0+800))]; % arbitrary point on line for creating vector
    d = abs(((p2(:,1)-p1(:,1)).*(p1(:,2)-p0(:,2))-(p2(:,2)-p1(:,2)).*(p1(:,1)-p0(:,1))) ./ sqrt((p2(:,1)-p1(:,1)).^2+(p2(:,2)-p1(:,2)).^2));
%     d_max = max(d); % identify maximum distance
    k_maxleft = find(ismembertol(d(1:(nan_loc-2)),d_max,0.000001)==1)'; % index of pt with max distance, numerical tolerance, nan_loc - 1 is same as 1
    k_maxright = find(ismembertol(d((nan_loc+2):end),d_max,0.000001)==1)'+nan_loc+1; % nan_loc+1 is same as end
    k_max = [k_maxleft; k_maxright]; % note: either the left or right k_max could be empty
    k_minleft = find(ismembertol(d(1:(nan_loc-2)),d_min,0.000001)==1)'; % index of pt with max distance, numerical tolerance
    k_minright = find(ismembertol(d((nan_loc+2):end),d_min,0.000001)==1)'+nan_loc+1;
    k_min = [k_minleft; k_minright]; % note: either the left or right k_min could be empty
    if length(k_max)>=2 % NEW 07/25/23 changed from ==2 to >= 2, for degenerate case where extra point within tolerance
        % eliminate point
        k_max_1 = k_max(1); % max d location for left poly
        k_max_2 = k_max(end); % max d location for right poly % NEW 07/25/23 changed from k_max(2) to k_max(end), for degenerate case where extra point within tolerance
%         Rch_cx = BBcx;
%         Rch_cy = BBcy;
%         Rch_cx = [BBcx1(1:k_min1); BBcx2(k_min2:end); BBcx1(end)];
%         Rch_cx(k_max_1+1:k_max_2-1)=[]; % remove points that are not part of convex hull, just for this specific region
%         Rch_cy(k_max_1+1:k_max_2-1)=[]; 
    end
    if length(k_max)==1 && ~ismembertol(BBcx(k_max),mask_x,0.000001)
        % add x_mask intercept as polygon point
        % k_xmask = find(ismembertol(BBcx(1:(nan_loc-2))),mask_x,0.000001)==1)';
        k_xmask = find(ismembertol(BBcx,mask_x,0.000001)==1)';
        y_max = 0;
        for ii=k_xmask
            if BBcy(ii)>y_max
               y_max = BBcy(ii);
            end
        end
        BBcx_alt = BBcx; BBcy_alt = BBcy;
        k_maxmask = k_xmask(find(BBcy(k_xmask)==y_max))'; 
        for ii=k_maxmask
           BBcy_alt(ii)=(BBcy(k_max)-(A(2,1)/A(2,2))*(mask_x-BBcx(k_max))); % where c3/c4 line intersects
        end
%         if length(k_maxmask)==1 % i.e., therefore the max y point and first entry for ccw rotation is from u_min
        k_max_1 = k_maxmask(1); % max distance location on left poly %%%%%%%%%don't want this to = 1 and end of left poly
        % k_max_1 = k_maxmask(end); % max distance location on left poly
        k_max_2 = k_max; % max distance location on right poly
%         else
%             k_max_2 = length(BBcx_alt); % max location should now be at beginning and end
%             k_max_1 = k_max(1);
%         end
%         Rch_cx = BBcx_alt;
%         Rch_cy = BBcy_alt;
%         Rch_cx(k_max_1+1:k_max_2-1)=[]; % remove points that are not part of convex hull, just for this specific region
%         Rch_cy(k_max_1+1:k_max_2-1)=[]; 
    end
    if length(k_max)==1 && ismembertol(BBcx(k_max),mask_x,0.000001)
        %remove all points b/w max x2 on mask and min x2 on mask (moving CCW)
        k_xmask = find(ismembertol(BBcx,mask_x,0.000001)==1)';
        y_max = 0;
        y_min = 15000;
        for ii=k_xmask
            if BBcy(ii)>=y_max % >= so that k_maxmask = end of right ploy not 1st entry
               y_max = BBcy(ii);
               k_maxmask = ii; 
            end
            if BBcy(ii)<y_min
               y_min = BBcy(ii);
               k_minmask = ii; 
            end
        end
        k_max_1 = k_minmask(end); % max distance location on left poly
        k_max_2 = k_maxmask; % max distance location on right poly
%         Rch_cx = BBcx;
%         Rch_cy = BBcy;
%         Rch_cx(k_maxmask+1:k_minmask-1)=[]; % remove points that are not part of convex hull, just for this specific region
%         Rch_cy(k_maxmask+1:k_minmask-1)=[];
    end
    % Same treatment for bottom of polygon
%     p0 = [Rch_cx, Rch_cy];
%     d = abs(((p2(:,1)-p1(:,1)).*(p1(:,2)-p0(:,2))-(p2(:,2)-p1(:,2)).*(p1(:,1)-p0(:,1))) ./ sqrt((p2(:,1)-p1(:,1)).^2+(p2(:,2)-p1(:,2)).^2));
%     d_min = min(d); % identify minimum distance
%     k_min = find(ismembertol(d(2:end),d_min,0.000001)==1)+1; % index of pt with max distance, numerical tolerance
    if length(k_min)>=2 % NEW 07/25/23 changed from ==2 to >= 2, for degenerate case where extra point within tolerance
        % eliminate point
        k_min_1 = k_min(1); % mim d location for left poly
        k_min_2 = k_min(end); % min d location for right poly % NEW 07/25/23 changed from k_min(2) to k_min(end), for degenerate case where extra point within tolerance
%         Rch_cx(k_min_1+1:k_min_2-1)=[]; % remove points that are not part of convex hull, just for this specific region
%         Rch_cy(k_min_1+1:k_min_2-1)=[]; 
    end
    if length(k_min)==1 && ~ismembertol(BBcy(k_min),mask_y,0.000001)
        % add y_mask intercept as polygon point
        k_ymask = find(ismembertol(BBcy,mask_y,0.000001)==1)';
        x_min = 5000;
        for ii=k_ymask
            if BBcx(ii)<x_min
               x_min = BBcx(ii);
            end
        end
        BBcx_alt = BBcx; BBcy_alt = BBcy;
        k_minmask = k_ymask(find(BBcx(k_ymask)==x_min))';
        for ii=k_minmask
           BBcx_alt(ii)=(BBcx(k_min)-(A(2,2)/A(2,1))*(mask_y-BBcy(k_min))); % where c3/c4 line intersects
        end
        k_min_1 = k_min; % minimum distance location on left poly
        k_min_2 = k_minmask(end); % minimum distance locatio on right poly
%         k_min_2 = k_minmask;
%         k_min_1 = k_min(1); % max location should now be at beginning and end
%         Rch_cx = Rchcx_alt;
%         Rch_cy = Rchcy_alt;
%         Rch_cx(k_min_1+1:k_min_2-1)=[]; % remove points that are not part of convex hull, just for this specific region
%         Rch_cy(k_min_1+1:k_min_2-1)=[]; 
    end
    if length(k_min)==1 && ismembertol(BBcy(k_min),mask_y,0.000001)
        %remove all points b/w max x2 on mask and min x2 on mask (moving CCW)
        k_ymask = find(ismembertol(BBcy,mask_y,0.000001)==1)';
%         if (k_ymask(end)-k_ymask(1))>1
%             x_max = 0;
%             x_min = 5000;
%             for ii=k_ymask
%                 if BBcx(ii)<x_min
%                    x_min = BBcx(ii);
%                    k_minmask = ii; %??
%                 end
%                 if BBcx(ii)>x_max
%                    x_max = BBcx(ii);
%                    k_maxmask = ii; %??
%                 end
%             end
%             Rch_cx = BBcx;
%             Rch_cy = BBcy;
%             Rch_cx(k_minmask+1:k_maxmask-1)=[]; % remove points that are not part of convex hull, just for this specific region
%             Rch_cy(k_minmask+1:k_maxmask-1)=[]; 
%         end
        k_min_1 = k_ymask(1); % minimum distance location on left poly
        k_min_2 = k_ymask(end); % minimum distance locatio on right poly
    end 
    %%%%%debug begin
    if length(k_min)>=3
        k_minleft = find(ismembertol(d(1:(nan_loc-2)),d_min,0.000001)==1)'; % index of pt with max distance, numerical tolerance
        k_minright = find(ismembertol(d((nan_loc+2):end),d_min,0.000001)==1)'+nan_loc+1;
    end
    if length(k_min)<=0
        k_minleft = find(ismembertol(d(1:(nan_loc-2)),d_min,0.000001)==1)'; % index of pt with max distance, numerical tolerance
        k_minright = find(ismembertol(d((nan_loc+2):end),d_min,0.000001)==1)'+nan_loc+1;
    end
    %%%%%% debug end
    Rch_cx = [BBcx(k_max_1:k_min_1); BBcx(k_min_2:k_max_2); BBcx(k_max_1)]; % left side of left poly and right side of right poly
    Rch_cy = [BBcy(k_max_1:k_min_1); BBcy(k_min_2:k_max_2); BBcy(k_max_1)];
end
