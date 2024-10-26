%% Task 4
%clear; clc; 

% Initialise parameter regimes
param_regimes = [1, 1/14, 1/20;
                 0.5, 1/14, 1/10;
                 1, 1/14, 1/200];

% Initialise Bridges
BsList = 10:10:100;

% Initialise how many combinations we want to run (testing purposes)
number_of_regimes = 3;
number_of_bridgetypes = 10;
simLength = 401;

% Run Simulation (testing - not using all bridges)
% [S, E, I, R] = spatialSim(param_regime_2, 10);

% Initialise Vectors
S = zeros(number_of_bridgetypes*number_of_regimes, simLength);
E = zeros(number_of_bridgetypes*number_of_regimes, simLength);
I = zeros(number_of_bridgetypes*number_of_regimes, simLength);
R = zeros(number_of_bridgetypes*number_of_regimes, simLength);

% Run Simulation
for i = 1:number_of_regimes
    for j = 1:number_of_bridgetypes
        k = j+(i-1)*10;
        [S(k,:), E(k,:), I(k,:), R(k,:)] = spatialSim(param_regimes(i,:), BsList(j));
    end
end

% Remove large number of windows
close all

%% 
% Plot Figures
for i = 1:number_of_regimes
    figure(i)
    clf
    for j = 1:number_of_bridgetypes
        k = j+(i-1)*10;
        subplot(2,5,j)
        hold on
        plot(0:1:(length(S)-1), S(k,:),'Color',[0.4660 0.6740 0.1880], 'LineWidth',1);
        plot(0:1:(length(S)-1), E(k,:),'Color',[0.9290 0.6940 0.1250], 'LineWidth',1);
        plot(0:1:(length(S)-1), I(k,:),'Color',[0.8500 0.3250 0.0980], 'LineWidth',1);
        plot(0:1:(length(S)-1), R(k,:),'Color',[0 0.4470 0.7410], 'LineWidth',1);
        xlabel("Time"); ylabel("Number of Agents")
        title(num2str(BsList(j))+ " Bridges")
        hold off
    end
    sgtitle("\alpha = "+num2str(param_regimes(i,1)) + ...
            ", \beta = "+num2str(round(param_regimes(i,2),3)) + ...
            ", \rho = "+num2str(param_regimes(i,3)))
    % Place legend in last plot
    legend('Susceptible', 'Exposed', 'Infectious', 'Recovered', 'Location','best')
    legend('boxoff')
end
