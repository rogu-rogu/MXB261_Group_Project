clear all
%Initial conditions 
S0 = 990; E0 = 0; I0 = 10; R0 = 0;
beta = 1/14;% Given
T = 360;  % Period
n_samples = 10;  % Number of LHS samples

lhs_samples_new = [rand(n_samples,1),rand(n_samples,1)];
alphas_new = lhs_samples_new(:, 1);  % Alpha samples
rhos_new = lhs_samples_new(:, 2);  % Rho samples

% Store successful parameter pairs
successful_params = [];
conditions = [];

for i = 1:n_samples
    % Extract sampled parameters
    alpha_new = alphas_new(i);
    rho_new = rhos_new(i);
    
    % Initialize initial state [S, E, I, R]
    X0 = [S0, E0, I0, R0];
    
    % Run SSA simulation 
    [t, X] = gillespieSSA(alpha_new, beta, rho_new, X0, T);
    
    % Check conditions
    S_T = X(1, end);  % Final susceptible count at time T
    if S_T >= S0 - 20 && S_T <= S0 - 10
        successful_params = [successful_params; alpha_new, rho_new];
        conditions = [conditions; 1];  % Condition (1)
    elseif S_T < S0 - 100
        successful_params = [successful_params; alpha_new, rho_new];
        conditions = [conditions; 2];  % Condition (2)
    end
end

figure;
hold on;
scatter(successful_params(conditions == 1, 1), successful_params(conditions == 1, 2), ...
        'r', 'DisplayName', 'Condition (1)');
scatter(successful_params(conditions == 2, 1), successful_params(conditions == 2, 2), ...
        'b', 'DisplayName', 'Condition (2)');
xlabel('\alpha');
ylabel('\rho');
legend('show');
title('Successful Parameter Combinations: \alpha vs \rho');
hold off;

%% Plotting Latin Hypercube
% % Only really works for low values of n_samples but hey! it makes sense!
% figure
% hold on 
% scatter(alphas_new, rhos_new,'.')
% grid on
% xticks(0:1/n_samples:1), xticklabels({})
% yticks(0:1/n_samples:1), yticklabels({})