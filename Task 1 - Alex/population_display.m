function population_display(t,X,method,alpha,beta,rho,xlimits,disp_type)
%{
inputs:
t - time range
X - matrix of population [ S(t) | E(t) | I(t) | R(t) ]
method - string with method that will be added to title
alpha, beta, rho - rate paramters for SEIR ODE
xlimits - upper limit for x axis
disp_type - show each sub-population on a separate sub-figure (0) or all on one
figure (1)
Output:
will produce the desired figure
%}
if nargin == 6
    usinglims = 0;
else
    usinglims = 1;
end
if nargin < 8
    disp_type = 2;
end
if disp_type == 1
    figure
    subplot(2,2,1)
    plot(t,X(:,1))
    title('subpopulation of Susceptible')
    subplot(2,2,2)
    plot(t,X(:,2))
    title('subpopulation of Exposed')
    subplot(2,2,3)
    plot(t,X(:,3))
    title('subpopulation of Infectious')
    subplot(2,2,4)
    plot(t,X(:,4))
    title('subpopulation of Recovered')
    ylabel('Subpopulation (people)')
    xlabel('Time (days)')
    sgtxt = [method ' method'];
    sgtitle(sgtxt)
    subtitle(" \alpha = " + alpha + ", \beta = " + beta + ", \rho = " + rho)
elseif disp_type == 2
    figure
    plot(t,X,'LineWidth',1)
    ylabel('Subpopulation (people)')
    xlabel('Time (days)')
    if usinglims
        xlim([0 xlimits])
    end
    legend('S','E','I','R')
    sgtxt = [method ' method'];
    title(sgtxt)
    subtitle(" \alpha = " + alpha + ", \beta = " + beta + ", \rho = " + rho)
end

end