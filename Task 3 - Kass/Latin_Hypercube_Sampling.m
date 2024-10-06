%% 3D Latin Hypercube Sampling
n = 16; % Number of samples
scale = (0:(1/n):1)'; scale(end) = []; 
%^^ Marking the leftmost value of n many segments ∈[0,1), the leftmost...
% value of Latin intervals
lhs_order = [randperm(n)', randperm(n)', randperm(n)']; 
%^^ Assigns all intervals one random number (unique to their hyperplane)

u = rand(n,3)/n; % Determines 3*n random numbers in a (1/n) square
parameters = u + scale(lhs_order); % Shifts values into their intervals
alpha_LHS = parameters(:,1);  
beta_LHS = parameters(:,2);  
rho_LHS = parameters(:,3);

%% 3D Orthogonal Sampling
sub_interval = 2; % Number of subspaces along a single dimension
subspace_count = sub_interval^3; % Number of subspaces in three dimensions
if mod(n,subspace_count) ~= 0
    error("The number of samples n must be divisible amongst "... 
       + string(subspace_count) + " subspaces (sub_interval^3).")
end
m = n / subspace_count; % Number of samples in a subspace

u = rand(m,3,subspace_count)/m; % Determines 3*m random numbers in...
% subspace_count many (1/m) squares


%% Plotting Latin Hypercube
% figure
% hold on 
% grid on
% scatter3(alpha, beta, rho)