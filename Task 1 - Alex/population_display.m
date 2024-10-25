function population_display(t,X,method,alpha,beta,rho,disp_type)
if nargin == 6
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
    legend('S','E','I','R')
    sgtxt = [method ' method'];
    title(sgtxt)
    subtitle(" \alpha = " + alpha + ", \beta = " + beta + ", \rho = " + rho)
end

end