clear
%% 3D Latin Hypercube Sampling
n = 27*2; % Number of samples

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
subinterval = 3; % Number of subspaces along a single dimension
subspace_count = subinterval^3; % Number of subspaces in three dimensions
if mod(n,subspace_count) ~= 0
    error("The number of samples n must be divisible amongst "... 
       + string(subspace_count) + " subspaces (sub_interval^3).")
end
m = n / subspace_count; % Number of samples in a subspace

u = rand(m,3,subspace_count); % Generating 3*m random numbers in...
%...subspace_count many squares ∈[0,1]  (in total 3*n random numbers)

% Distributing samples to Latin Hypercubes--------------------------------]
internal_scale = (0:1/m:1)'; internal_scale(end) = []; 
internal_order = zeros(m,3);
subspace = u/m; % Shrink to ∈[0, 1/m]
for i = 1:subspace_count
    for dimension = 1:3 % Assign intervals per LHS
        internal_order(:,dimension) = randperm(m)';
    end
    subspace(:,:,i) = subspace(:,:,i) + internal_scale(internal_order);
    %^^ Shifts values into their intervals
end

% Distributing subspaces to Latin Hypercubes------------------------------]
subspace_scale = (0:1/subinterval:1)'; subspace_scale(end) = []; 
subspace_order = [randperm(subinterval)', randperm(subinterval)',...
    + randperm(subinterval)']; % Only performed over one space, so no loop
space = subspace/subinterval; % Shrink to ∈[0, 1/subinterval]
parameters_O = zeros(n,3);
for i = 1:subspace_count
    % Converting from (m x 3 x subspace_count) to (m*subspace_count x 3)
    space2list = [m*(i-1)+1, m*i];
end

%% Plotting
% Only really works for low values of n_samples but hey! it makes sense!
figure
hold on 
scatter3(alpha_LHS, beta_LHS, rho_LHS, '.')
grid on
xticks(0:1/n:1), xticklabels({})
yticks(0:1/n:1), yticklabels({})
