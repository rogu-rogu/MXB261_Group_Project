%% 3D Latin Hypercube Sampling
n = 10; % Number of samples
u = rand(n,3)/n; % Determining 3*n random numbers in a (1/n) square
scale = (0:(1/n):1)'; scale(end) = []; 
%^^ Marking the leftmost value of n many segments ∈[0,1), the leftmost...
% value of Latin intervals
lhs_order = [randperm(n)', randperm(n)', randperm(n)']; 
%^^ Assigns all intervals one random number (unique to their hyperplane)
parameters = u + scale(lhs_order); % Shifts values into their intervals
alpha = parameters(:,1);  beta = parameters(:,2);  rho = parameters(:,3);

%% Orthogonal Sampling


%% Plotting Latin Hypercube
% figure
% hold on 
% grid on
% scatter3(alpha, beta, rho)