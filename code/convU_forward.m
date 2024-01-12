function [Rch_cx, Rch_cy] = convU_forward(BBpgon, A, mask_x, mask_y)
% Remove points between "tips" of u_max and u_min regions from polyshape because intermediated controls can reach from that space (this is taking convex hull of just that one region)
% Note: this should really apply at bottom tips of reach sets too, butnot worth coding because those aren't realistic part of state space and won't even appear in figures
[BBcx,BBcy] = boundary(BBpgon);
[BBcx,BBcy] = poly2ccw(BBcx,BBcy); %counterclockwise
BBcx(abs(BBcx)<0.1)=0;
BBcy(abs(BBcy)<0.1)=0;
maxy = max(BBcy); % identify maximum y coordinate - the tip
miny = min(BBcy);
k_max = find(BBcy==maxy)';
minx_maxy = 5000;
maxx_maxy = 0;

% Find points with max and min y value (to connect), with updated indices
if (sum(isnan(BBcx))==0) % only one polygon present
    for ii=k_max
        if (BBcx(ii)<minx_maxy)
            minx_maxy=BBcx(ii);
            k_topleft = ii;
        end
    end
    if ~(k_topleft==1) % that is, does not begin and end with maximum y point
        BBcx = [BBcx(k_topleft:end); BBcx(2:k_topleft)]; % change so that do begin and end with maximum y point
        BBcy = [BBcy(k_topleft:end); BBcy(2:k_topleft)];
    end
    k_max = find(BBcy==maxy)';
    k_min = find(BBcy==miny)';
    Rch_cx = BBcx;
    Rch_cy = BBcy;
    Rch_cx(k_max(2)+1:k_max(3)-1)=[]; % connect top points by removing pts between
    Rch_cy(k_max(2)+1:k_max(3)-1)=[]; 
    Rch_cx(k_min(1)+1:k_min(2)-1)=[]; % connect bottom points by removing pts between - 
    % based on ordering (above) these indices are before those for "top" so don't have to re-index
    Rch_cy(k_min(1)+1:k_min(2)-1)=[]; 
else
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
    if k_topleft<nan_loc
        if ~(k_topleft==1) % that is, does not begin and end with maximum y point
            BBcx1 = [BBcx(k_topleft:(nan_loc-1)); BBcx(2:k_topleft)]; % change so that do begin and end with maximum y point
            BBcy1 = [BBcy(k_topleft:(nan_loc-1)); BBcy(2:k_topleft)];
        end
        if ~(k_topright==(nan_loc+1)) % that is, does not begin and end with maximum y point
            BBcx2 = [BBcx(k_topright:end); BBcx((nan_loc+2):k_topright)]; % change so that do begin and end with maximum y point
            BBcy2 = [BBcy(k_topright:end); BBcy((nan_loc+2):k_topright)];
        end
    else
        if ~(k_topleft==1) % that is, does not begin and end with maximum y point
            BBcx1 = [BBcx(k_topleft:end); BBcx((nan_loc+2):k_topleft)]; % change so that do begin and end with maximum y point
            BBcy1 = [BBcy(k_topleft:end); BBcy((nan_loc+2):k_topleft)];
        end
        if ~(k_topright==(nan_loc+1)) % that is, does not begin and end with maximum y point
            BBcx2 = [BBcx(k_topright:(nan_loc-1)); BBcx(2:k_topright)]; % change so that do begin and end with maximum y point
            BBcy2 = [BBcy(k_topright:(nan_loc-1)); BBcy(2:k_topright)];
        end  
    end
    k_min1 = find(BBcy1==miny)';
    k_min2 = find(BBcy2==miny)';
    Rch_cx = [BBcx1(1:k_min1); BBcx2(k_min2:end); BBcx1(end)]; % left side of left poly and right side of right poly
    Rch_cy = [BBcy1(1:k_min1); BBcy2(k_min2:end); BBcy1(end)];
end