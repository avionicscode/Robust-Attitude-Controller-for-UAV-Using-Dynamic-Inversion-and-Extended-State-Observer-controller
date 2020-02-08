t=0:10/(length(t)-1):10;
T=0:10/(length(T)-1):10;
subplot(3,1,1)
plot(T,y(:,1),'r-.') % phi
hold on
plot(t,x(:,1),'b') % phi
title('Angle de Roulis')
ylabel('\phi °','fontsize',22)
xlabel('temps s')
legend('sans ESO','avec ESO')
subplot(3,1,2)
plot(T,y(:,4),'r-.')  %teta
hold on
plot(t,x(:,4),'b') % teta
title('Angle de tangage')
ylabel('\theta °','fontsize',22)
xlabel('temps s')
legend('sans ESO','avec ESO')
subplot(3,1,3)
plot(T,y(:,7),'r-.')  %psi
hold on
plot(t,x(:,7),'b') % psi
title('Angle de lacet')
ylabel('\psi °','fontsize',22)
xlabel('temps s')
legend('sans ESO','avec ESO')