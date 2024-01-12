alpha_1 = 0.12; % equilibrium ratio of internal P to dissolved+adsorbed P  [-] %0.12
alpha_2 = 0.004; % fraction of dissolved+adsorbed pool that is dissolved at equilibrium [-] % 0.004
alpha_3 = 0.365; % rate of internal diffusion[1/T] % 0.365
alpha_4 = 45; % maximum plant uptake of dissolved P [P] % 30
alpha_5 = 40; % Michaelis-Menten parameter for plant uptake of dissolved P [P] % 15
alpha_6 = 0.1; % fraction of streambed P lost by erosion per unit flow [-] ??????? % 0.05
beta_1 = 0.0025;  % fraction of dissolved P lost with runoff per unit flow [-] % 0.008
beta_2 = 0.00025;  % fraction of particulate P lost with runoff per unit flow [-] (for 1 m depth, Q_bar=4, 0.00025 -> 1.0 mm loss) % 0.001
beta_3 = 0.35;  % fraction of fixed P in water column that settles into streambed [-] ????????
Q_bar = 4; % mean flow ?????????

x_vec = [218; 1500; 1.5]; % X(1,:) = [ 280; 1400; 4] or [ 120; 700; 5]?

A_s = zeros(3,3); % linear, stochastic processes
A_s(1,1) = -(beta_1*alpha_2 + beta_2*(1 - alpha_2));
A_s(2,2) = -beta_2;
A_s(3,1) = 0;%alpha_3*(beta_1*alpha_2 + beta_2*(1 - alpha_2));%got rid of flow effect for dilution (offsets)
A_s(3,2) = beta_3*beta_2;
A_s(3,3) = -alpha_6;

Sigma_vec=A_s*x_vec;
NonLinSigma = zeros(3,3);
NonLinSigma(1,1) = Sigma_vec(1)^2;
NonLinSigma(2,2) = Sigma_vec(2)^2;
NonLinSigma(3,3) = Sigma_vec(3)^2;
NonLinSigma(1,2) = 0.95*sqrt(NonLinSigma(1,1)*NonLinSigma(2,2));%*(sign(Sigma_vec(1))*sign(Sigma_vec(2)));
NonLinSigma(2,1) = NonLinSigma(1,2);
NonLinSigma(1,3) = 0.95*sqrt(NonLinSigma(1,1)*NonLinSigma(3,3));%*(sign(Sigma_vec(1))*sign(Sigma_vec(3)));
NonLinSigma(3,1) = NonLinSigma(1,3);
NonLinSigma(2,3) = 0.95*sqrt(NonLinSigma(2,2)*NonLinSigma(3,3));%*(sign(Sigma_vec(2))*sign(Sigma_vec(3)));
NonLinSigma(3,2) = NonLinSigma(2,3);

mu = [0 0 0];
Sigma = NonLinSigma;
% Sigma = [1 2.99; 2.99 9];

x1 = -1:0.05:1;
x2 = -1:0.05:1;
x3 = -1:0.05:1;
[X1,X2, X3] = meshgrid(x1,x2,x3);
X = [X1(:) X2(:) X3(:)];

y = mvnpdf(X,mu,Sigma);
y = reshape(y,length(x3),length(x2),length(x1));

figure
surf(x1,x2,y(:,:,21))
caxis([min(y(:))-0.5*range(y(:)),max(y(:))])
axis([-3 3 -3 3 0 0.4])
xlabel('x1')
ylabel('x2')
zlabel('Probability Density')

% c = Sigma_vec(1)+Sigma_vec(2)+Sigma_vec(3);
% d = Sigma_vec(1)*Sigma_vec(2)*Sigma_vec(3);
% sig_roots = roots([2 -c 0 d]);
% sig_roots = sort(sig_roots(imag(sig_roots)==0));
% num_roots = size(sig_roots);
% max_D = sqrt(Sigma_vec(1)*Sigma_vec(2));
% if num_roots >= 2
%     if max_D >= sig_roots(2)
%         off_D = sig_roots(2) - 0.01;
%     else
%         off_D = max_D - 0.01;
%     end
% else
%     off_D = max_D - 0.01;
% end

