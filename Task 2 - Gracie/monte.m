clear all
%Using Task 1 Code for ODE comparison. 
function dXdt = Epidemic_ode(t,X,a,b,r)
dSdt = -a*X(1)*X(3)/(X(1)+X(2)+X(3)+X(4));
dEdt = a*X(1)*X(3)/(X(1)+X(2)+X(3)+X(4)) - b*X(2);
dIdt = b*X(2) - r*X(3);
dRdt = r*X(3);

dXdt = [dSdt,dEdt,dIdt,dRdt]';
end
function [t_span, X_euler] = Euler_method_Epidemic(t_span, X_IC, alpha, beta, rho, h)
  %range of time for simulation
X_euler = zeros(length(t_span),4); %initialise
X_euler(1,:) = X_IC;

for i = 1:length(t_span)-1
    dXdt = Epidemic_ode(t_span,X_euler(i,:),alpha,beta,rho);
    X_euler(i+1,:) = X_euler(i,:) + h*dXdt';
end

end
%Monte Carlo Function 
function monte_carlo_simulation(Nsim)
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
    t_span = 0:0.5:360;
    X_IC = [990,0,10,0];
    [t_euler, X_euler] = Euler_method_Epidemic(t_span,X_IC,alpha,beta,rho,0.5);
       
    

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
            X(:,ii+1) = max(X(:,ii) + nu*(Y'), 0);
            t(ii+1) = t(ii) + dt;
        end
        D{k} = [t;X];
    end
    toc
    
    X_sum = zeros(N, Nt);
    
    for k = 1:Nsim
        X_sum = X_sum + D{k}(2:end, :) ;
         
    end
    % Calculate expected value
    X_est = X_sum / Nsim;

    figure;
    % Plot Susceptible
    plot(t, X_est(1, :), 'lineWidth', 2, 'Color', [0 0.4470 0.7410]);  % Solid blue line for Susceptible
    hold on;
    plot(t_euler, X_euler(:, 1), '--', 'lineWidth', 2, 'Color', [0 0.4470 0.7410]);  % Dashed blue line for Susceptible ODE
    
    % Plot Exposed
    plot(t, X_est(2, :), 'lineWidth', 2, 'Color', [0.8500 0.3250 0.0980]);  % Solid red line for Exposed
    plot(t_span, X_euler(:, 2), '--', 'lineWidth', 2, 'Color', [0.8500 0.3250 0.0980]);  % Dashed red line for Exposed ODE
    
    % Plot Infected
    plot(t, X_est(3, :), 'lineWidth', 2, 'Color', [0.8 0.6 0]);  % Solid dark yellow line for Infected
    plot(t_euler, X_euler(:, 3), '--', 'lineWidth', 2, 'Color', [0.8 0.6 0]);  % Dashed dark yellow line for Infected ODE
    
    % Plot Recovered
    plot(t, X_est(4, :), 'lineWidth', 2, 'Color', 'm');  % Solid magenta line for Recovered
    plot(t_euler, X_euler(:, 4), '--', 'lineWidth', 2, 'Color', 'm');  % Dashed magenta line for Recovered ODE
    
    xlabel('Time (days)');
    ylabel('Population');
    legend('Susceptible', 'Susceptible ODE', 'Exposed', 'Exposed ODE', 'Infected', 'Infected ODE', 'Recovered', 'Recovered ODE','Location', 'southoutside', 'NumColumns', 2);

    title(['Monte Carlo Integration of SEIR Model using Tau-Leap with Nsim = ', num2str(Nsim)]);
    
    hold off;

end
%Produce the plots 
monte_carlo_simulation(3)
monte_carlo_simulation(6)
monte_carlo_simulation(12)
monte_carlo_simulation(25)
monte_carlo_simulation(50)
monte_carlo_simulation(100)