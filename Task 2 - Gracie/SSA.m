clear all
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

% propensity functions
a = @(X,t) [alpha*X(1)*X(3)/sum(X),beta*X(2),rho*X(3)];
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

% ODE steady state
%need code from part 1 c
t_span = 0:0.5:360;
X_IC = [990,0,10,0];
[t_euler, X_euler] = Euler_method_Epidemic(t_span,X_IC,alpha,beta,rho,0.5);
%Ss;
%Is;
%Es;
%Rs;


%% simulate with SSA and plot
%% Code adapated from week 8 interative lecture content. 
Nsim = 3;

T = 360;
D = cell(Nsim,1);
b = 1;
tic;
for k=1:Nsim
    % initialise t and X
    t = zeros(1,b);
    X = zeros(N,b);
    X(:,1) = X0;
    ii = 1;
    while true
        % update propensities
        a_tmp = a(X(:,ii),t(ii));
        a0 = sum(a_tmp);
        % sample the waiting time dt ~ Exp(a0)
        u1 = rand;
        dt = -log(u1)/a0;
        if t(end) + dt > T
            % stop SSA and plot output
            break;
        else
            % sample next reaction with lookup method
            P = cumsum(a_tmp./a0);
            u2 = rand;
            i=1;
            while u2 > P(i)
                i = i+1;
            end
            j = i;
        end
        ii = ii + 1;
        if ii > length(t)
            t = [t,zeros(1,b)];
            X = [X,zeros(N,b)];
        end
        % update state and time
        t(1,ii) = t(1,ii-1) + dt; 
        X(:,ii) = X(:,ii-1) + nu(:,j);
        %disp(['time = ',num2str(t(1,ii))]);
    end
    D{k} = [t;X];
end
toc
%% plot the simulation results
figure();
sgtitle("Gillespie stochastic simulation algorithm (SSA) comparison to ODE")
 for k=1:Nsim
 t = D{k}(1,:);
 X = D{k}(2:end,:);
             subplot(2,2,1);
             hold on;
             plot(t,X(1,:),'lineWidth',1);
             ylabel('S(t)');
             xlabel('t')
             
             subplot(2,2,2);
             hold on;
             plot(t,X(2,:),'lineWidth',1);
             xlabel('t');
             ylabel('E(t)');
             subplot(2,2,3)
             hold on
             plot(t,X(3,:),'lineWidth',1);
             xlabel('t');
             ylabel('I(t)');
             subplot(2,2,4)
             hold on
             plot(t,X(4,:),'lineWidth',1);
             xlabel('t');
             ylabel('R(t)');
 end
 for k=1:Nsim
             subplot(2,2,1);
             hold on;
             plot(t_span,X_euler(:,1).*ones(size(t_span)),'--k','lineWidth',2);
             subplot(2,2,2);
             hold on;
             plot(t_span,X_euler( :,2).*ones(size(t_span)),'--k','lineWidth',2);
             subplot(2,2,3);
             hold on;
             plot(t_span,X_euler(:, 3).*ones(size(t_span)),'--k','lineWidth',2);
             subplot(2,2,4)
             hold on
             plot(t_span,X_euler(:,4).*ones(size(t_span)),'--k','lineWidth',2);

 end

 %%
 %%