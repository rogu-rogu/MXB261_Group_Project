%% Orthogonal Sampling
% Informed by the procedure developed in "Populations of models, 
% Experimental Designs and coverage of parameter space by Latin Hypercube 
% and Orthogonal Sampling"; K. Burrage, P. Burrage, D. Donovan, B. Thompson
% doi: 10.1016/j.procs.2015.05.383

% Initialising -----------------------------------------------------------]
d = 3; % Number of dimensions
n = 27; % Number of samples, number of subspaces
p = nthroot(n,d); % Number of subspace intervals
if floor(p)~=p, error("n must have a positive integer cube root."), end
m = n/p; % Number of sample intervals (within subspaces)
subscale = 0:(1/p):1; subscale(end) = []; % Marks the beginning position...
%...of p many Latin intervals in ∈[0,1]. Positions subspaces in space
samplescale = 0:(1/m):1; samplescale(end) = []; % As above, but marks...
%...m many Latin intervals. Positions samples within a subspace

% Constructs an organised list of subspaces ------------------------------]
subco = 1:p; % Range of subspace positions within the space
subspaces = zeros(n,d); % List of subspace indexes
subspaces(:,1) = repelem(subco,m);              % Changes every m values
subspaces(:,2) = repmat(repelem(subco,p),1,p);  % Changes every p values
subspaces(:,3) = repmat(subco,1,m);             % Changes every value

% Isolates hyperplanes of subspaces --------------------------------------]
hyperplanes = false(n,d,p);
for i = subco, hyperplanes(:,:,i) = subspaces == i; end
% Locates i in subspaces, isolating the hyperplane where parameter == i

% Assigns samples within subspaces according to LHS ----------------------]
samplespaces = zeros(n,d);
for i = subco % Loops through hyperplanes of subpaces
    for j = 1:d % Loops through parameters
        val = randperm(m); % Assign random positions with no duplicates
        samplespaces(hyperplanes(:,j,i),j) = val;
    end
end

% Generates random samples -----------------------------------------------]
u = rand(n,d); % Generates uniform random numbers (3*n for 3D)
u_subspace = u/m; % Shrinks values into a (1/m) hypercube
u_subspace = u_subspace + samplescale(samplespaces); 
% Shifts values into their intervals, positioning samples within a subspace
u_space = u_subspace/p; % Shrinks values into a subspace ((1/m) hypercube)
u_space = u_space + subscale(subspaces); % Shifts subspaces into position
% Generates outputs ------------------------------------------------------]
alpha = u_space(:,1);
beta = u_space(:,2);
rho = u_space(:,3);

figure
hold on 
scatter3(alpha, beta, rho, '.')
grid on
xticks(0:1/n:1), xticklabels({})
yticks(0:1/n:1), yticklabels({})
view(3)