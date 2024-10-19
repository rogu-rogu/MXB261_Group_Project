
function monte_carlo_simulation(Nsim)
    % Monte Carlo integration using modified Tau-Leap code

   %Monte Carlo integration using modified Tau Leap code from week 11
    %interactive lecture 
    % initial state
    X0 = [990,0, 10, 0];

    % parameters
    alpha = 1/5;
    beta = 1/14;
    rho = 1/10;
    M = 3;
    N = 4;
    % stochiometric matrix
    nu = [-1 0 0  ;
          1 -1 0 ;
          0 1 -1;
           0 0 1];
    % population
    P = 1000;

    % propensity functions
    a = @(X,t) [alpha*X(1)*X(3)/sum(X),beta*X(2),rho*X(3)];

    
    T = 360;   
    Nt = 721;  
    dt = T / (Nt - 1); 

    
    D = cell(Nsim, 1);

    % Perform simulations
    %% simulate with tau-leap and plot
    
    
    T = 360;
    Nt = 721;
    D = cell(Nsim,1);
    
    tic;
    
    for k=1:Nsim
        X = zeros(N,Nt);
        t = zeros(1,Nt);
        dt = T/(Nt-1);
        X(:,1) = X0;
        t(1) = 0;
        for ii=1:(Nt-1)
            % evaluate propensities
            a_tmp = a(X(:,ii),t(ii));
            % generate M Poisson random variables
            Y = poissrnd(a_tmp*dt);
            % update state and time
            X(:,ii+1) = X(:,ii) + nu*(Y'); 
            t(ii+1) = t(ii) + dt;
        end
        D{k} = [t;X];
    end
    toc

    
    X_sum = zeros(N, Nt);

    for k = 1:Nsim
        X_sum = X_sum + D{k}(2:end, :); 
    end

    % Calculate expected value
    X_est = X_sum / Nsim;

    figure
    % Plot 
    plot(t, X_est(1, :), 'lineWidth', 2);  % E Susceptible
    hold on
    plot(t, X_est(2, :), 'lineWidth', 2); % E Exposed
    hold on
    plot(t, X_est(3, :), 'lineWidth', 2);    % E Infected
    hold on
    plot(t, X_est(4, :), 'lineWidth', 2);  % E Recovered
    hold on

    xlabel('Time (days)');
    ylabel('Population');
    legend('Susceptible', 'Exposed', 'Infected', 'Recovered');
    title(['Monte Carlo Simulation of SEIR Model using Tau-Leap with Nsim = ', num2str(Nsim)]);
    hold off;

end
monte_carlo_simulation(3)
monte_carlo_simulation(6)
monte_carlo_simulation(12)
monte_carlo_simulation(25)
monte_carlo_simulation(50)
monte_carlo_simulation(100)