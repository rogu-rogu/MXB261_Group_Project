clear,clc
alpha = 1/5;
beta = 1/14;
rho = 1/10;
t_span = 0:0.5:360;
X_IC = [990,0,10,0];
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
population_display(t_euler,X_euler,'Euler')
population_display(t_ode,X_ode,'ode45')
X_diff = X_ode - X_euler;
population_display(t_ode,X_diff,'difference between')
%% range of alpha 
alpha = 0.025:(0.05-0.025):1;
%will be simulating using ode45 as a guess.
S = zeros(1,length(alpha));
for i = 1:length(alpha)
   [t_ode,X_ode] = ode45(@(t,X) Epidemic_ode(t,X,alpha(i),beta,rho),t_span,X_IC);
   S(i) = X_ode(end,1);
end
plot(alpha,S)
title('relationship between S and alpha')
subtitle('after 360 days')
xlabel('alpha value')
ylabel('Susceptible subpopulation')