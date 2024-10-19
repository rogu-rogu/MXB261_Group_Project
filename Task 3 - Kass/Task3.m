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
partial_LHS = sample_success_LHS > 0; % Has at least one successful trial

%% Implement Orthogonal Sampling on parameter space (alpha, beta, rho)
orthog = Orthogonal_Sampling(n,d);
sample_success_orthog = condition_test(orthog, m, X0, T);
success_orthog = sample_success_orthog >= 1 - tol;
partial_orthog = sample_success_orthog > 0;

%% Generating figures
c = flip(hot); % Colour gradient of successful trials 
c(1:90,:) = []; % Removing lightest colours, hard to see against white
gr = '#0bcd08'; % Successful sample colour
gr2 = '#099ec2'; % Alternate successful sample colour


% FIGURE 1: Plot test results of LHS samples =============================]
f1 = figure; colormap(c), hold on, grid on, axis([0 1 0 1 0 1]), axis equal
ax1 = gca;
%----Plots-------------------------------------%
LHS_o = LHS(success_LHS,:); % Successful samples
LHS_x = LHS(~success_LHS,:); % Unsuccessful samples
scatter3(LHS_o(:,1), LHS_o(:,2), LHS_o(:,3), 'o', ...
    'MarkerEdgeColor', gr);
scatter3(LHS_x(:,1), LHS_x(:,2), LHS_x(:,3), [], ...
    sample_success_LHS(~success_LHS), 'x', 'LineWidth', 1);
%----Labels---------------------%
title("Latin Hypercube Sampling")
xlabel("alpha"), ylabel("beta"), zlabel("rho")
lgd = legend("Successful samples", "Failed samples " + newline + ...
    "(darker colour" + newline + "corresponds to" + newline + ...
    "greater proportion" + newline + "of successful trials)");
%----Positioning------------%
lgd.Location = "eastoutside";
%=========================================================================]


% FIGURE 2: Plot test results of orthogonal samples ======================]
f2 = figure; colormap(c), hold on, grid on, axis([0 1 0 1 0 1]), axis equal
ax2 = gca;
%----Plots-------------------------%
orthog_o = orthog(success_orthog,:);
orthog_x = orthog(~success_orthog,:);
scatter3(orthog_o(:,1), orthog_o(:,2), orthog_o(:,3), ...
    'MarkerEdgeColor', gr);
scatter3(orthog_x(:,1), orthog_x(:,2), orthog_x(:,3), [], ...
    sample_success_orthog(~success_orthog), 'x', 'LineWidth', 1)
%----Labels----------------%
title("Orthogonal Sampling")
xlabel("alpha"), ylabel("beta"), zlabel("rho")
lgd = legend("Successful samples", "Failed samples " + newline +...
    "(darker colour" + newline + "corresponds to" + newline +...
    "greater proportion" + newline + "of successful trials)");
%----Positioning------------%
lgd.Location = "eastoutside";
%=========================================================================]


% FIGURE 3: Plot test results of Latin Hypercube and orthogonal samples ==]
f3 = figure; colormap(c), hold on, grid on, axis([0 1 0 1 0 1]), axis equal
ax3 = gca;
%----Plots---------------------------------------------------%
scatter3(LHS_o(:,1), LHS_o(:,2), LHS_o(:,3), ... % LHS success
    'MarkerEdgeColor', gr); 
scatter3(orthog_o(:,1), orthog_o(:,2), orthog_o(:,3), ... % orthog success
    'MarkerEdgeColor', gr2);
scatter3(LHS_x(:,1), LHS_x(:,2), LHS_x(:,3), [], ... % LHS fail
    sample_success_LHS(~success_LHS), 'x', 'LineWidth', 1);
scatter3(orthog_x(:,1), orthog_x(:,2), orthog_x(:,3), [], ... % orthog fail
    sample_success_orthog(~success_orthog), 'x', 'LineWidth', 1)
%----Labels------------------------------------%
title("Latin Hypercube and Orthogonal Sampling")
xlabel("alpha"), ylabel("beta"), zlabel("rho")
lgd = legend("LHS successful samples", "Orthogonal successful " + ...
    newline + "samples", "Failed samples (darker" + newline +  ...
    "colour corresponds to" + newline + "greater proportion" + newline +...
    "of successful trials)");
%----Positioning------------%
lgd.Location = "eastoutside";
%=========================================================================]


% FIGURE 4: Convex hull of successful LHS samples ========================]
f4 = figure; colormap(c), hold on, grid on, axis([0 1 0 1 0 1]), axis equal
ax4 = gca;
%----Plots------------------------------------------------------------%
LHS_box_success = convhulln(LHS_o); % Box containing successful samples
LHS_p = LHS(partial_LHS,:); % Samples with at least one successful trial
LHS_box_partial = convhulln(LHS_p); %^^^Box containing samples as above
trisurf(LHS_box_success, LHS_o(:,1), LHS_o(:,2), LHS_o(:,3), ...
    'FaceColor', gr);
trisurf(LHS_box_partial, LHS_p(:,1), LHS_p(:,2), LHS_p(:,3), ...
    'FaceColor', gr2, 'FaceAlpha', 0.2)
%----Labels--------------------------------------%
title("Parameter space of successful LHS samples")
xlabel("alpha"), ylabel("beta"), zlabel("rho")
lgd = legend("Polygon containing" + newline + "successful samples", ...
    "Polygon containing" + newline + "samples with at" + newline + ...
    "least one" + newline + "successful trial");
%----Positioning------------%
lgd.Location = "eastoutside";
%=========================================================================]


% Figure 5: Convex hull of successful orthogonal samples =================]
f5 = figure; colormap(c), hold on, grid on, axis([0 1 0 1 0 1]), axis equal
ax5 = gca;
%----Plots------------------------------%
orthog_box_success = convhulln(orthog_o);
orthog_p = orthog(partial_orthog,:);
orthog_box_partial = convhulln(orthog_p);
trisurf(orthog_box_success, orthog_o(:,1), orthog_o(:,2), orthog_o(:,3),...
    'FaceColor', gr);
trisurf(orthog_box_partial, orthog_p(:,1), orthog_p(:,2), orthog_p(:,3),...
    'FaceColor', gr2, 'FaceAlpha', 0.2)
%----Labels--------------------------------------------%
title("Parameter space of successful orthognal samples")
xlabel("alpha"), ylabel("beta"), zlabel("rho")
lgd = legend("Polygon containing" + newline + "successful samples", ...
    "Polygon containing" + newline + "samples with at" + newline + ...
    "least one" + newline + "successful trial");
%----Positioning------------%
lgd.Location = "eastoutside";
%=========================================================================]

%% Visualisation
% Default 3-dimensional view
view(ax1, 3), saveas(f1, '3D_view/Figure_1.png')
view(ax2, 3), saveas(f2, '3D_view/Figure_2.png')
view(ax3, 3), saveas(f3, '3D_view/Figure_3.png')
view(ax4, 3), saveas(f4, '3D_view/Figure_4.png')
view(ax5, 3), saveas(f5, '3D_view/Figure_5.png')

% 2D view of alpha and beta
view(ax1, 0, 90), saveas(f1, 'alpha_beta/Figure_1.png')
view(ax2, 0, 90), saveas(f2, 'alpha_beta/Figure_2.png')
view(ax3, 0, 90), saveas(f3, 'alpha_beta/Figure_3.png')
view(ax4, 0, 90), saveas(f4, 'alpha_beta/Figure_4.png')
view(ax5, 0, 90), saveas(f5, 'alpha_beta/Figure_5.png')

% 2D view of alpha and rho
view(ax1, 0, 0), saveas(f1, 'alpha_rho/Figure_1.png')
view(ax2, 0, 0), saveas(f2, 'alpha_rho/Figure_2.png')
view(ax3, 0, 0), saveas(f3, 'alpha_rho/Figure_3.png')
view(ax4, 0, 0), saveas(f4, 'alpha_rho/Figure_4.png')
view(ax5, 0, 0), saveas(f5, 'alpha_rho/Figure_5.png')

% 2D view of beta and rho
view(ax1, 90, 0), saveas(f1, 'beta_rho/Figure_1.png')
view(ax2, 90, 0), saveas(f2, 'beta_rho/Figure_2.png')
view(ax3, 90, 0), saveas(f3, 'beta_rho/Figure_3.png')
view(ax4, 90, 0), saveas(f4, 'beta_rho/Figure_4.png')
view(ax5, 90, 0), saveas(f5, 'beta_rho/Figure_5.png')