%% Create safe set
% 97.5% chance of satisfying water quality, satisfy nutrient deficit constraint
% Def_lim = 0.75; % nutrient deficit limit (i.e. 75% of max)
% Load_lim = 0.6*Q_bar*((beta_1*alpha_2+beta_2*(1-alpha_2))*x1_init + beta_2*x2_init); % 226 and 1850 are estimated IC
x3_nom = 0.6; %0.5; 
ps_nom = 0.0; %0.4;
x1_min = Def_lim/(1-Def_lim)*alpha_5;
x2_max = (Load_lim - Q_bar*(beta_1*alpha_2+beta_2*(1-alpha_2))*x1_min) / (Q_bar*beta_2);
x1_max = Load_lim /(Q_bar*(beta_1*alpha_2+beta_2*(1-alpha_2)));
cx = [x1_min, x1_min, x1_max]; % corner points, x-coord
cy = [0, x2_max, 0]; % corner points, y-coord
poly_points = [cx; cy]';
pgon=polyshape(poly_points);
SafeSet=pgon;

%% Create equilibrium line
x1_eq = [0 250]; % 120-130 is relevant portion for base
x2_eq = alpha_3*x1_eq/(beta_2*Q_bar+alpha_3*alpha_1);
% figure
% hold on
% plot(SafeSet,'FaceColor',[0 0 0], 'FaceAlpha',1);
% plot(x1_eq,x2_eq,'r');


%% Check for intersection
[x1i,x2i] = polyxpoly(x1_eq,x2_eq,cx,cy);


