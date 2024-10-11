%% Initialising
% Variables set within this file are stored here for easy access
X0 = [990  0  10  0]; % Initial simulation conditions
T = 360; % Maximum time span
d = 3;   % Number of parameters/dimensions (alpha, beta, rho)
n = 5^d; % Number of samples
m = 10;  % Number of trials
tol = 0; % Tolerance for proportion of failed trials, (success >= 1 - tol)

%% Implement Latin Hypercube Sampling on parameter space (alpha, beta, rho)
LHS = Latin_Hypercube_Sampling(n,d);
sample_success_LHS = condition_test(LHS, m, X0, T);
success_LHS = sample_success_LHS >= 1 - tol;

%% Implement Orthogonal Sampling on parameter space (alpha, beta, rho)
orthog = Orthogonal_Sampling(n,d);
sample_success_orthog = condition_test(orthog, m, X0, T);
success_orthog = sample_success_orthog >= 1 - tol;

%% Plotting
c = flip(hot); % Colour gradient of successful trials 
c(1:50,:) = []; % Removing lightest colours, hard to see against white


% Figure 1: Plot test results of LHS samples------------------------------]
f1 = figure; colormap(c), hold on, grid on, axis([0 1 0 1 0 1]), axis equal
%----Plots-------------------------------------%
LHS_o = LHS(success_LHS,:); % Successful samples
LHS_x = LHS(~success_LHS,:); % Unsuccessful samples
scatter3(LHS_o(:,1), LHS_o(:,2), LHS_o(:,3), 'o', 'MarkerFaceColor', 'g');
scatter3(LHS_x(:,1), LHS_x(:,2), LHS_x(:,3), [], ...
    sample_success_LHS(~success_LHS), 'x', 'LineWidth', 1);
%----Labels---------------------%
title("Latin Hypercube Sampling")
xlabel("alpha"), ylabel("beta"), zlabel("rho")
lgd = legend("Successful samples", "Failed samples " + newline + ...
    "(darker colour" + newline + "corresponds to" + newline + ...
    "greater proportion" + newline + "of successful trials)");
%----Positioning---------------%
f.Position = [30, 80, 700, 500];
lgd.Location = "eastoutside";
view(3) %-----------------------------------------------------------------]

% Figure 2: Plot test results of orthogonal samples-----------------------]
f2 = figure; colormap(c), hold on, grid on, axis([0 1 0 1 0 1]), axis equal
%----Plots-------------------------%
orthog_o = orthog(success_orthog,:);
orthog_x = orthog(~success_orthog,:);
scatter3(orthog_o(:,1), orthog_o(:,2), orthog_o(:,3), ...
    'MarkerFaceColor','g');
scatter3(orthog_x(:,1), orthog_x(:,2), orthog_x(:,3), [], ...
    sample_success_orthog(~success_orthog), 'x', 'LineWidth', 1)
%----Labels---------------------%
title("Orthogonal Sampling")
xlabel("alpha"), ylabel("beta"), zlabel("rho")
lgd = legend("Successful samples", "Failed samples " + newline +...
    "(darker colour" + newline + "corresponds to" + newline +...
    "greater proportion" + newline + "of successful trials)");
%----Positioning---------------%
f.Position = [30, 80, 700, 500];
lgd.Location = "eastoutside";
view(3) %-----------------------------------------------------------------]