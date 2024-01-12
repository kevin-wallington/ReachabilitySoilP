%% Initialize vectors and set controls
horizon=450;
X = cell(12,1);
O = cell(12,1);
X(:) = {zeros(horizon, 2)}; % state: dissolved+adsorbed, internally fixed, streambed
% Initial x1  = 52 g/Mg (Bray P ave) * 1.5 (convert to total ads+dis) * 1.45 Mg/m^3 * 10,000 m^2/ha = 1131 kg/ha (assumes modelling 1 m depth)(226 for .2 m)
% X{1}(1,:)=[ 250; 2000; 1.45]; X{2}(1,:)=[ 300; 1250; 1.45]; X{3}(1,:)=[ 500; 790; 1.45]; X{4}(1,:)=[ 700; 500; 1.45]; X{5}(1,:)=[ 50; 1750; 1.45]; % for mining, u=u_min
% x1_init = 226;
X{1}(1,:)=[ 300;1800]; X{2}(1,:)=[ 300; 1200]; X{3}(1,:)=[ 300; 840]; X{4}(1,:)=[ 300; 480]; X{5}(1,:)=[ 30; 1750]; % for mining, u=u_min
X{6}(1,:)=[ 30; 1500]; X{7}(1,:)=[ 30; 1280]; X{8}(1,:)=[ 30; 973]; X{9}(1,:)=[ 30; 680]; X{10}(1,:)=[ 150; 200]; % for increase, u=u_max
X{11}(1,:)=[ x1_init; x2_init]; X{12}(1,:)=[ x1_init; x2_init]; % estimated IC
O(:) = {zeros(horizon, 2)}; % output: adsorbed pool, river export

%% Simulate and plot sample trajectory
for kk=1:12
    if(kk<=5)
        u = u_min; % mining/drawing, natural weathering only
    elseif(kk<=10)
        u = u_max; % build-up, 30 taken as max fertilizer
    end
    if(kk==11); u = 21; end % status quo % +7 if including subsoil supple
    if(kk==12); u = u_min; end % optimal traj from estimate IC (change to max U later inside loop)
    for ii=2:horizon
        if(kk==12)
            if(ii==40); u = u_max; end %switch to build up, max fert old-40 , new-44
            if(ii>=43); u= 16; end % stabilize inside target  old-43 , new-46 % +7 if including subsoil supply
        end
        f_nonlin = [(-alpha_4*(X{kk}(ii-1,1)/(X{kk}(ii-1,1)+alpha_5))); 0];
        b = [u;0]; 
        X{kk}(ii,:) = (A * X{kk}(ii-1,:)' + f_nonlin + b)'; % update state
        O{kk}(ii-1,:) = C * X{kk}(ii-1,:)';    
    end
end

%% Plot state space evolution
figure
hold on
plot(SafeSet,'FaceColor',[0 0 0], 'FaceAlpha',1);
plot(SafeInv,'FaceColor',[0.5 0.5 0.5], 'FaceAlpha',1); 
for kk=2:5
    plot(X{kk}(:,1), X{kk}(:,2), 'color', "#A2142F",'linewidth', 2) %"#A2142F" for mining/drawdown, add '-x' for points
    plot([X{kk}(1,1)], [X{kk}(1,2)],'o','MarkerFaceColor','black','MarkerEdgeColor','black', 'MarkerSize',10,'linewidth',3);
end
for kk=7:10
    plot(X{kk}(:,1), X{kk}(:,2), 'color', "#EDB120",'linewidth', 2) %"#EDB120" for build up/increase
    plot([X{kk}(1,1)], [X{kk}(1,2)],'o','MarkerFaceColor','black','MarkerEdgeColor','black', 'MarkerSize',10,'linewidth',3);
end
    plot([X{10}(1,1)], [X{10}(1,2)],'o','MarkerFaceColor','white','MarkerEdgeColor','white', 'MarkerSize',10,'linewidth',3);
% for kk=11
%     plot(X{kk}(:,1), X{kk}(:,2), 'color', "#0072BD",'linewidth', 2) %"#0072BD" for status quo
% end
% for kk=12
%     plot(X{kk}(1:43,1), X{kk}(1:43,2), 'color', "#77AC30", 'linewidth', 3) %"#77AC30" for optimal trajectory, old-43 , new-46
% end
for kk=2
    plot(X{kk}(1:23,1), X{kk}(1:23,2), 'color', "#77AC30", 'linewidth', 4) %"#77AC30" for bang-bang trajectory, old-43 , new-46
end
for kk=8
    plot(X{kk}(2:3,1), X{kk}(2:3,2), 'color', "#77AC30", 'linewidth', 4) %"#77AC30" for bang-bang trajectory, old-43 , new-46
end
plot([X{2}(1,1)], [X{2}(1,2)],'o','MarkerFaceColor','black','MarkerEdgeColor','black', 'MarkerSize',10,'linewidth',3);
% plot formatting
xlim([0 305])
ylim([0 1800])
ax = gca; % axes handle
ax.FontSize = 18; 
title('', 'fontsize', 26)
xlabel('Labile P (kg/ha)', 'fontsize', 26)
ylabel('Stable P (kg/ha)', 'fontsize', 26)
% legend('', 'Location', 'Northeast')
% legend boxoff

% plot initial point
% plot([x1_init],[x2_init],'pentagram','MarkerFaceColor','black','MarkerEdgeColor','black', 'MarkerSize',50,'linewidth',2);
