%% Initialize vectors and set controls
horizon = 5000;
X = zeros(horizon, 2); % state: dissolved+adsorbed, internally fixed, streambed
% Initial x1  = 52 g/Mg (Bray P ave) * 1.5 (convert to total ads+dis) * 1.45 Mg/m^3 * 10,000 m^2/ha = 1131 kg/ha (assumes modelling 1 m depth)(226 for .2 m)
X(1,:)=[ x1_init; x2_init]; % estimated IC
O = zeros(horizon, 2); % output: adsorbed pool, river export
hit_time = 0;

%% Target and switching region
pgon=SafeInv; % Safe_Inv created by SimpleModel_Invariance.m
polyout = polybuffer(pgon,0.1);
[cx,cy] = boundary(polyout); 

%% Simulate and plot sample trajectory
fert_rate = 21 + u_min; % fertilizer rate [kg/ha] 
for ii=2:horizon
    f_nonlin = [(-alpha_4*(X(ii-1,1)/(X(ii-1,1)+alpha_5))); 0];
    b = [fert_rate;0]; 
    X(ii,:) = (A * X(ii-1,:)' + f_nonlin + b)'; % update state
    O(ii-1,:) = C * X(ii-1,:)';  
    if(inpolygon(X(ii,1),X(ii,2), cx,cy) && hit_time == 0) % in target region for first time
        hit_time = ii-1;
    end
end