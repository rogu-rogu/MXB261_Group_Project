function dXdt = Epidemic_ode(t,X,a,b,r)
dSdt = -a*X(1)*X(3)/(X(1)+X(2)+X(3)+X(4));
dEdt = a*X(1)*X(3)/(X(1)+X(2)+X(3)+X(4)) - b*X(2);
dIdt = b*X(2) - r*X(3);
dRdt = r*X(3);

dXdt = [dSdt,dEdt,dIdt,dRdt]';
end