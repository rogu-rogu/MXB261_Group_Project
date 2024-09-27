function [t_span, X_euler] = Euler_method_Epidemic(t_span, X_IC, alpha, beta, rho, h)
  %range of time for simulation
X_euler = zeros(length(t_span),4); %initialise
X_euler(1,:) = X_IC;

for i = 1:length(t_span)-1
    dXdt = Epidemic_ode(t_span,X_euler(i,:),alpha,beta,rho);
    X_euler(i+1,:) = X_euler(i,:) + h*dXdt';
end

end