%% Initialize vectors and arrays
x1_invlimit = 1000;
pgon = SafeInv;
[cx,cy] = boundary(pgon);
x2_val = max(cy); % evaluate at this x2 value

%% Back-propagate
u= u_min;
k = 1;
for jj=1:100 % x1 value
        c5 = A(1,1) - A(1,2) * A(2,1) / A(2,2);
        c6 = (A(1,2) / A(2,2)) * (x2_val) + u;
        a = c5; 
        b = c5*alpha_5 - alpha_4 + c6 - (jj);
        c = c6*alpha_5 - alpha_5*(jj);
        root_out = roots([a b c]);
        if (isreal(root_out))
            x1_invlimit = min(x1_invlimit, jj);
        end
end
mask_x = x1_invlimit; % testing *0.5
%% Plot
% figure
% scatter(x1_sp,x2_sp,2,back_x1);
% %caxis([0, max_minx1]);
% xlim([0 100])
% ylim([0 1500])
% ax = gca; % axes handle
% ax.FontSize = 18; 
% title('', 'fontsize', 26)
% xlabel('Short-term P (kg/ha)', 'fontsize', 26)
% ylabel('Long-term P (kg/ha)', 'fontsize', 26)
% a = colorbar;
% a.Label.String = 'Minimum solution to inverse map';