function [t, X] = gillespieSSA(alpha, beta, rho, X0, T)
% GILLESPIESSA - Simulate a stochastic SEIR model using Gillespie algorithm
% The SEIR model assumes a well-mixed, fixed volume P with four distinct
% species (compartments: S = "susceptible", E = "exposed", I =
% "infectious", R = "recovered") and a set of three reactions (exposure, 
% infection, recovery).
%
%   [t, X] = gillespieSSA(alpha, beta, rho, X0, T) simulates an epidemic
%   model with initial conditions X0 in domain t∈[0,T]. The parameters
%   alpha, beta and rho affect the rate of reaction, or state change
%   probability.
%
%   Input Arguments
%     alpha - scalar between 0 and 1, double
%       Parameter affecting the propensity function of exposure; the rate
%       of change from S to E
%     beta - scalar between 0 and 1, double
%       Paramater affecting the propensity function of infection (referring
%       to the end of the latent period); the rate of change from E to I
%     rho - scalar between 0 and 1, double
%       Parameter affecting the propensity function of recovery; the rate
%       of change from I to R.
%     X0 - (1x4) uint16 vector <---- actually double but we could change this to marginally save memory
%       The initial condition, X(t=0) = [S,E,I,R].
%     T - scalar, double
%       Determines the maximum extent of time the simulation computes,
%       modelling from t=0 to t=T.
%
%   Output Arguments
%     t - (1x?) double matrix
%       Time values associated with each state X(t).
%     X - (4x?) double matrix
%       State vector over time (row 1 = S, row 2 = E, row 3 = I, row 4 = R)

    S = X0(1); E = X0(2); I = X0(3); R = X0(4); % Initial conditions 
    X = X0'; 
    t = 0;
    nu = [-1 0 0; 1 -1 0; 0 1 -1; 0 0 1];  % Stoichiometric vectors

    while t(end) < T
        P = S + E + I + R;  % Total population
        a = [alpha * S * I / P, beta * E, rho * I];  % Propensity functions
        a0 = sum(a);
        
        % Sample the waiting time
        dt = -log(rand) / a0;
        new_t = t(end) + dt; 
        if new_t > T; break;
        else t(end+1) = new_t; end
        
        % Determine which reaction occurs
        p = cumsum(a) / a0;
        u = rand; % u ~ Uniform(0,1)
        j = find(u < p, 1);
        
        % Update state
        X = [X, X(:,end) + nu(:,j)];
        S = X(1, end); E = X(2, end); I = X(3, end); R = X(4, end);
    end
end