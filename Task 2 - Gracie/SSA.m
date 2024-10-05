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

% ODE steady state
%need code from part 1 c
%Ss;
%Is;
%Es;
%Rs;


%% simulate with SSA and plot
%% Code adapated from week 8 interative lecture content. 
Nsim = 10;

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
 for k=1:Nsim
 t = D{k}(1,:);
 X = D{k}(2:end,:);
            subplot(2,2,1);
             hold on;
             plot(t,X(1,:),'lineWidth',1);
             ylabel('S');
             xlabel('time')
             
             subplot(2,2,2);
             hold on;
             plot(t,X(2,:),'lineWidth',1);
             xlabel('time');
             ylabel('E');
             subplot(2,2,3)
             hold on
             plot(t,X(3,:),'lineWidth',1);
             xlabel('time');
             ylabel('I');
             subplot(2,2,4)
             hold on
             plot(t,X(4,:),'lineWidth',1);
             xlabel('time');
             ylabel('R');
 end
%need code from part 1 c
 %for k=1:Nsim
 %subplot(1,2,1);
             %hold on;
             %plot(t,As*ones(size(t)),'--k','lineWidth',2);
             
             %subplot(1,2,2);
             %hold on;
             %plot(t,Bs*ones(size(t)),'--k','lineWidth',2);
 %end
 %subplot(1,2,1);
 %set(gca, 'YScale', 'log')