function [S, E, I, R] = spatialSim(param_regime, B)

% Initialise Parameters
alpha = param_regime(1);
beta = param_regime(2);
rho = param_regime(3);

% Initialise Location
gridxLength = 101; 
gridyLength = 101;
riverLocation = 51;

% Initialise Population parameters
S0 = 990; 
E0 = 0; 
I0 = 10; 
R0 = 0;
P = S0 + I0 + E0 + R0;

% Set maximum simulation Time
maxSimTime = 100;

% Initialise state trackers
S = zeros(1, maxSimTime);
E = zeros(1, maxSimTime);
I = zeros(1, maxSimTime);
R = zeros(1, maxSimTime);

% Create grid
gridSpace = zeros(gridxLength, gridyLength);
gridSpace(riverLocation,:) = 1;

% Initialise agent position and state
posStateMatrix = zeros(P,3);
posStateMatrix(:,3) = 1;

% Populate Grid
agentCount = 0;
I0Set = 0;
while agentCount < P
    x = randi(gridxLength);
    y = randi(gridyLength);
    if gridSpace(y,x) == 0
        agentCount = agentCount+1;
        gridSpace(y,x) = 1;
        posStateMatrix(agentCount, [1,2]) = [x, y];
        if all(I0Set < I0 & y > riverLocation)
            I0Set = I0Set+1;
            posStateMatrix(agentCount, 3) = 3;
        end
    end
end

% Add bridges
bridgeLocations = randperm(gridxLength);
gridSpace(riverLocation, bridgeLocations(1:B)) = 0;

% define time step
deltat = 1;

% Setup video
writerObj = VideoWriter(num2str(alpha)+"regime"+num2str(B)+"bridges.avi");
writerObj.FrameRate = 10;
open(writerObj);

figure('Position', [50 50 900 700])
hold on

for i = 1:maxSimTime
    % Move agents
    for j = 1:P
        direction = rand;
        if direction < 0.25 % Move North
            if (posStateMatrix(j,2)-1) >= 1
                if gridSpace(posStateMatrix(j,2)-1, posStateMatrix(j,1)) == 0
                    gridSpace(posStateMatrix(j,2)-1, posStateMatrix(j,1)) = 1;
                    gridSpace(posStateMatrix(j,2), posStateMatrix(j,1)) = 0;
                    posStateMatrix(j,2) = posStateMatrix(j,2)-1;
                end
            end
        elseif direction < 0.5 % Move South
            if (posStateMatrix(j,2)+1) <= 101
                if gridSpace(posStateMatrix(j,2)+1, posStateMatrix(j,1)) == 0
                    gridSpace(posStateMatrix(j,2)+1, posStateMatrix(j,1)) = 1;
                    gridSpace(posStateMatrix(j,2), posStateMatrix(j,1)) = 0;
                    posStateMatrix(j,2) = posStateMatrix(j,2)+1;
                end
            end
        elseif direction < 0.75 % Move East
            if (posStateMatrix(j,1)+1) <= 101
                if gridSpace(posStateMatrix(j,2), posStateMatrix(j,1)+1) == 0
                    gridSpace(posStateMatrix(j,2), posStateMatrix(j,1)+1) = 1;
                    gridSpace(posStateMatrix(j,2), posStateMatrix(j,1)) = 0;
                    posStateMatrix(j,1) = posStateMatrix(j,1)+1;
                end
            end
        else % Move West
            if (posStateMatrix(j,1)-1) >= 1
                if gridSpace(posStateMatrix(j,2), posStateMatrix(j,1)-1) == 0
                    gridSpace(posStateMatrix(j,2), posStateMatrix(j,1)-1) = 1;
                    gridSpace(posStateMatrix(j,2), posStateMatrix(j,1)) = 0;
                    posStateMatrix(j,1) = posStateMatrix(j,1)-1;
                end
            end
        end
    end
    
    % Set agent State
    oldStates = posStateMatrix(:,3);
    for j = 1:P
        if oldStates(j) == 2 % For exposed 
            % convert to Infected
            if rand <= beta
                posStateMatrix(j,3) = 3;
            end
        elseif oldStates(j) == 3 % For infected
            %  Infect a neighbour
            if rand <= alpha
                infectDirection = rand;
                if infectDirection < 0.25 % Infect North
                    findPosState = [posStateMatrix(j,1) posStateMatrix(j, 2)-1 1];
                    findPosStateIndex = find(ismember(posStateMatrix,findPosState,'rows'));
                elseif infectDirection < 0.5 % Infect South
                    findPosState = [posStateMatrix(j,1) posStateMatrix(j, 2)+1 1];
                    findPosStateIndex = find(ismember(posStateMatrix,findPosState,'rows'));
                elseif infectDirection < 0.75 % infect East
                    findPosState = [posStateMatrix(j,1) posStateMatrix(j, 2)+1 1];
                    findPosStateIndex = find(ismember(posStateMatrix,findPosState,'rows'));
                else % Infect West
                    findPosState = [posStateMatrix(j,1) posStateMatrix(j, 2)+1 1];
                    findPosStateIndex = find(ismember(posStateMatrix,findPosState,'rows'));
                end
                if findPosStateIndex
                    posStateMatrix(findPosStateIndex, 3) = 2;
                end
            end
            % Recover
            if rand <= rho
                posStateMatrix(j,3) = 4;
            end
        end
    end

    % Record state counts
    S(i) = nnz(posStateMatrix(:,3) == 1);
    E(i) = nnz(posStateMatrix(:,3) == 2);
    I(i) = nnz(posStateMatrix(:,3) == 3);
    R(i) = nnz(posStateMatrix(:,3) == 4);

    % Plot video
    clf
    hold on
    plot(posStateMatrix(posStateMatrix(:,3) == 1, 1), posStateMatrix(posStateMatrix(:,3) == 1, 2), ...
        'o','LineWidth',4,'MarkerSize',4,'Color',[0.4660 0.6740 0.1880])
    plot(posStateMatrix(posStateMatrix(:,3) == 2, 1), posStateMatrix(posStateMatrix(:,3) == 2, 2), ...
        'o','LineWidth',4,'MarkerSize',4,'Color',[0.9290 0.6940 0.1250])
    plot(posStateMatrix(posStateMatrix(:,3) == 3, 1), posStateMatrix(posStateMatrix(:,3) == 3, 2), ...
        'o','LineWidth',4,'MarkerSize',4,'Color',[0.8500 0.3250 0.0980])
    plot(posStateMatrix(posStateMatrix(:,3) == 4, 1), posStateMatrix(posStateMatrix(:,3) == 4, 2), ...
        'o','LineWidth',4,'MarkerSize',4,'Color',[0 0.4470 0.7410])
    legend('Susceptible', 'Exposed', 'Infectious', 'Recovered')
    xlabel("X Position"); ylabel("Y Position")
    xlim([0 gridxLength]) % Set grid size
    ylim([0 gridyLength])
    yregion(riverLocation-0.5, riverLocation+0.5)
    box on
    frame = getframe;
    writeVideo(writerObj,frame);

end

% end video
close(writerObj)

% add initial values to trackers
S = [S0 S];
E = [E0 E];
I = [I0 I];
R = [R0 R];

% Plot populations
hold on
plot(0:deltat:maxSimTime, S);
plot(0:deltat:maxSimTime, E);
plot(0:deltat:maxSimTime, I);
plot(0:deltat:maxSimTime, R);
xlabel("Time"); ylabel("Number of Agents")
title("\alpha = "+num2str(alpha)+", \beta = "+num2str(beta)+", \rho = "+num2str(rho), ", "+num2str(B)+" Bridges")
legend('Susceptible', 'Exposed', 'Infectious', 'Recovered')

end