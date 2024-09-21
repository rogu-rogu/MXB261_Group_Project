%% Task 4
% rename this using source control! the option is right beneath 
% "Extract Conflict Markers to File" and above "Move"

% Initialise parameter regimes
param_regime_1 = [1,1/14,1/20];
param_regime_2 = [0.5,1/14,1/10];

% Initialise Bridges
BsList = 10:10:100;

% Run Simulation (testing - not using all bridges)
[S, E, I, R] = spatialSim(param_regime_2, 10);


