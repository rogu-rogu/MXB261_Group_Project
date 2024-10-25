clear all
%Initial conditions 
S0 = 990; E0 = 0; I0 = 10; R0 = 0;
beta = 1/14;% Given
T = 360;  % Period
n_samples = 1000;  % Number of LHS samples

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
    
    % Initialise initial state [S, E, I, R] 
    X0 = [S0, E0, I0, R0];
    
    % Run SSA simulation 
    [t, X] = gillespieSSA_SEIR(alpha_new, beta, rho_new, X0, T);
    
    % Check conditions
    S_T = X(1, end);  % Final susceptible count at time T
    S_10 = X(1,10);
    if S_T >= S0 - 20 && S_T <= S_10 - 10
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

% Gillespie SSA function
function [t, x] = gillespieSSA_SEIR(alpha, beta, rho, X0, T)
    % Gillespie Stochastic Simulation Algorithm for an SEIR model
    % Inputs:
    %   alpha - infection rate parameter
    %   beta - incubation rate parameter
    %   rho - recovery rate parameter
    %   X0 - initial state [S0, E0, I0, R0]
    %   T - total simulation time
    % Outputs:
    %   t - time vector
    %   x - state vector with each row [S, E, I, R] over time

    %% Initialise variables
    t(1) = 0; % Initialise the simulation time
    count = 0; % Initialise count variable
    x(:, 1) = X0; % Initialise storage of abundances with initial state
    S = X0(1); E = X0(2); I = X0(3); R = X0(4); % Initial states 
    % Define reaction rate constants
   

    %% Begin SSA Loop
    while max(t) < T
        % Update count variable
        count = count + 1;

        % Calculate propensities based on reactions
        N = S + E + I + R; % Total population at current state
        a(1) = alpha * S * I / N; 
        a(2) = beta * E;         
        a(3) = rho * I;        
        
        % Total propensity
        a0 = sum(a);

        % Calculate the time step
        u1 = rand; 
        delta_t = -1 / a0 * log(u1); % Calculate time step
        t(count + 1) = t(count) + delta_t; % Update time

        % Determine which reaction occurs
        u2 = rand; % Another random number
        c = cumsum(a); % Calculate cumulative propensities
        if u2 * a0 <= c(1) 
            v = [-1; 1; 0; 0];
        elseif u2 * a0 <= c(2)
            v = [0; -1; 1; 0];
        elseif u2 * a0 <= c(3) 
            v = [0; 0; -1; 1];
        end

        % Update
        x(:, count + 1) = x(:, count) + v;

        % Update
        S = x(1, count + 1); E = x(2, count + 1); I = x(3, count + 1); R = x(4, count + 1);
    end
end