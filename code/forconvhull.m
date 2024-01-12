[cx,cy] = boundary(pgon); % will wrap back to original pt to close polygon
pgon_hull = convhull(cx,cy);
num_pts = length(pgon_hull);
hull_points = zeros(num_pts,2);
for kk=1:num_pts
    hull_points(kk,1)= cx(pgon_hull(kk));
    hull_points(kk,2)= cy(pgon_hull(kk));
end
pgon=polyshape(hull_points);