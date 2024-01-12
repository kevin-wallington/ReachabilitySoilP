x1_min = Def_lim/(1-Def_lim)*alpha_5;
x2_max = (Load_lim - Q_bar*(beta_1*alpha_2+beta_2*(1-alpha_2))*x1_min) / (Q_bar*beta_2);
x1_max = Load_lim /(Q_bar*(beta_1*alpha_2+beta_2*(1-alpha_2)));

A = [x1_min, x1_max];
B = [x2_max 0];

xlim = [0 200];
m = (B(2)-B(1))/(A(2)-A(1));
n = -A(2)*m + B(2);
y1 = m*xlim(1) + n;
y2 = m*xlim(2) + n;
hold on
line([xlim(1) xlim(2)],[y1 y2],'Color', 'yellow')

% And a line for the crop nutrient deficit
line([x1_min x1_min],[0 2000],'Color', 'red')
