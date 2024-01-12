%% Initialize vectors and arrays
x1_sp = zeros(1500000,1);
x2_sp = zeros(1500000,1);
back_x1 = zeros(1500000,1);
backprop_min = zeros(1500,1000);
max_minx1 = -1000;

%% Back-propagate
u= u_max;
k = 1;
for ii=1:1500 %linearize plant uptake at every boundary point then back propagate each point
    for jj=1:1000
            c5 = A(1,1) - A(1,2) * A(2,1) / A(2,2);
            c6 = (A(1,2) / A(2,2)) * (ii) + u;
            a = c5; 
            b = c5*alpha_5 - alpha_4 + c6 - (jj);
            c = c6*alpha_5 - alpha_5*(jj);
            root_out = roots([a b c]);
            if (isreal(root_out))
                backprop_min(ii,jj)= min(root_out);
                back_x1(k)= min(root_out);
                max_minx1 = max(max_minx1, min(root_out));
            else
                backprop_min(ii,jj)= NaN;
                back_x1(k)= NaN;
            end
            x1_sp(k) = jj;
            x2_sp(k) = ii;
            k=k+1;
    end
end

%% Plot
figure
scatter(x1_sp,x2_sp,2,back_x1);
%caxis([0, max_minx1]);
xlim([0 1000])
ylim([0 1500])
ax = gca; % axes handle
ax.FontSize = 18; 
title('', 'fontsize', 26)
xlabel('Short-term P (kg/ha)', 'fontsize', 26)
ylabel('Long-term P (kg/ha)', 'fontsize', 26)
a = colorbar;
a.Label.String = 'Minimum solution to inverse map';