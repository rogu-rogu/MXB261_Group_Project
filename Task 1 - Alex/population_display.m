function population_display(t,X,method)
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
sgtxt = [method ' method'];
sgtitle(sgtxt)

end