%% Orthogonal Sampling
% Informed by the procedure developed in "Populations of models, 
% Experimental Designs and coverage of parameter space by Latin Hypercube 
% and Orthogonal Sampling"; K. Burrage, P. Burrage, D. Donovan, B. Thompson
% doi: 10.1016/j.procs.2015.05.383

d = 3; % Number of dimensions
n = 27; % Number of samples, number of subspaces
p = nthroot(n,d); % Length of subspace intervals
m = n/p; % Length of sample intervals (within subspaces)
if floor(p)~=p, error("n must have an integer cube root."), end

% Constructs an organised list of subspaces ------------------------------]
subco = 1:p; % Range of subspace positions within the space
subspaces = zeros(n,d); % List of subspace indexes
subspaces(:,1) = repelem(subco,m);              % Changes every m values
subspaces(:,2) = repmat(repelem(subco,p),1,p);  % Changes every p values
subspaces(:,3) = repmat(subco,1,m);             % Changes every value

% Isolates hyperplanes to assign samples according to LHS ----------------]
hyperplanes = false(n,d,p);
for i = subco, hyperplanes(:,:,i) = subspaces == i; end
% Locates i in subspaces, isolating the hyperplane where parameter == i

% Randomly orders unique sample positions in each hyperplane -------------]
sampleco = 1:m; % Range of sample positions within subspaces


% 
% for i = subco 
%     set = subco == i;
%     x(:,1,i) = repelem(set,m); % Changes every m values
%     x(:,2,i) = repmat(repelem(set,p),1,p); % Changes every p values
%     x(:,3,i) = repmat(set,1,m); % Changes every value
% end

% Constructs the structure of the space: logical matrices...
    %...of the position of coordinates in an ordered arrangement. This
    % essentially is an ordered list of subspaces; which permutations can
    % then be applied to, to rearrange samples inside and to rearrange
    % their position in the space.