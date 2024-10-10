clear
%% 3D Latin Hypercube Sampling
n = 16; % Number of samples

% Demarcates the beginning of n many Latin intervals, in ∈[0,1]
scale = (0:(1/n):1); scale(end) = []; 
% Assigns all intervals one random number (unique to their hyperplane)
lhs_order = [randperm(n); randperm(n); randperm(n)]'; 

u = rand(n,3)/n; % Determines 3*n random numbers in a (1/n) cube
parameters_LHS = u + scale(lhs_order); % Shifts values into their intervals
alpha_LHS = parameters_LHS(:,1);  
beta_LHS = parameters_LHS(:,2);  
rho_LHS = parameters_LHS(:,3);


%% Plotting
% Only really works for low values of n_samples but hey! it makes sense!
figure
hold on 
scatter3(parameters(:,1), parameters(:,2), parameters(:,3), '.')
grid on
xticks(0:1/n:1), xticklabels({})
yticks(0:1/n:1), yticklabels({})
view(3)