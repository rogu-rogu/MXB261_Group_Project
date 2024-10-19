clear,clc
%% equilibrium analysis
syms S E I a P b r L
Xdashsyms = [(-a*S*I/P); (a*S*I/P - b*E); (b*E - r*I)]; %ODEs
Jsyms = jacobian(Xdashsyms, [S,E,I]); %jacobian
eigen_syms = -Jsyms + L*eye(3); %equation to find lambda values
eigen_0 = subs(eigen_syms,I,0);
lambda_sol = solve(det(eigen_0)==0,L) %get the three eigenvalues
%% initialise
alpha = 1/5;
beta = 1/14;
rho = 1/10;
t_span = 0:0.5:360;
X_IC = [990,0,10,0];
P = sum(X_IC);
%% Euler
[t_euler, X_euler] = Euler_method_Epidemic(t_span,X_IC,alpha,beta,rho,0.5);
population_display(t_euler,X_euler,'Euler')
%% Euler vs MATLABs ode45
X_IC = [990,0,10,0];
%ode45 expects y(t,y) for the input order of the function
[t_ode,X_ode] = ode45(@(t,X) Epidemic_ode(t,X,alpha,beta,rho),t_span,X_IC);
population_display(t_ode,X_ode,'ode45')
%% alpha = 10
alpha = 10;
X_IC = [990,0,10,0];
[t_euler, X_euler] = Euler_method_Epidemic(t_span,X_IC,alpha,beta,rho,0.5);
[t_ode,X_ode] = ode45(@(t,X) Epidemic_ode(t,X,alpha,beta,rho),t_span,X_IC);
%display Euler Solution
population_display(t_euler,X_euler,'Euler')
%display ode45 ssolution
population_display(t_ode,X_ode,'ode45')
%display X1 - X2, individual difference within S,E,I,R
population_display(t_ode,X_ode-X_euler,'Difference of each')
%convert into a vector
X_diff = vecnorm(X_ode - X_euler,1,2);
figure
plot(t_euler,X_diff)
ylabel('||X_{ode45} - X_{euler}||')
xlabel('Time (days)')
title('Norm of difference per day')
%% range of alpha 
alpha_range = 0.025:(0.05-0.025):1;
%will be simulating using ode45 as a guess.
S = zeros(1,length(alpha_range));
for i = 1:length(alpha_range)
   [t_ode,X_ode] = ode45(@(t,X) Epidemic_ode(t,X,alpha_range(i),beta,rho),t_span,X_IC);
   S(i) = X_ode(end,1);
end
plot(alpha_range,S)
title('relationship between S and alpha')
subtitle('after 360 days')
xlabel('alpha value')
ylabel('Susceptible subpopulation')
%% extra stuff that will be removed
rho_range = 0.025:(0.05-0.025):1;
alpha_range = 0.025:(0.05-0.025):1;
R = zeros(length(alpha_range),length(rho_range));
for j = 1:length(rho_range)
    for i = 1:length(alpha_range)
        [t_ode,X_ode] = ode45(@(t,X) Epidemic_ode(t,X,alpha_range(i),beta,rho_range(j)),t_span,X_IC);
        R(i,j) = X_ode(end,1);
    end
end
figure
mesh(alpha_range,rho_range,R)
xlabel('\rho')
ylabel('\alpha')
zlabel('S(360)')
title('Affect of transmission and recovery rate on final S sub-pop')