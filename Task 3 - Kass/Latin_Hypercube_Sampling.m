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

%% 3D Orthogonal Sampling
d = 3; % Number of dimensions
n = 27; % Number of samples, number of subspaces
m = nthroot(n,d); % Length of subinterval. Cube root for 3 dimensions
x = n/m; % Number of subspaces in hyperplane
if floor(m)~=m
    error("n must have an integer cube root.")
end

scale = (0:(1/n):1); scale(end) = []; 
% Demarcates the beginning of m many Latin intervals, in ∈[0,1)
subscale = (0:(1/m):1); subscale(end) = []; 
sampleorder = zeros(d,n);
suborder = zeros(d,n);

% for i = 1:x
%     focus = m*(i-1)+1 : m*i;
%     for parameter = 1:d
% 
%         sampleorder(parameter,focus) = randperm(m);
% 
%         if parameter == 1, val = 1:m;
%         elseif parameter == 2, val = mod(i-1,m)+1;
%         else, val = floor((i-1)/m)+1;
%         end
%         suborder(parameter,focus) = val;
%     end
% end



%% bug testing


%% Plotting
% Only really works for low values of n_samples but hey! it makes sense!
figure
hold on 
scatter3(parameters(:,1), parameters(:,2), parameters(:,3), '.')
grid on
xticks(0:1/n:1), xticklabels({})
yticks(0:1/n:1), yticklabels({})
view(3)