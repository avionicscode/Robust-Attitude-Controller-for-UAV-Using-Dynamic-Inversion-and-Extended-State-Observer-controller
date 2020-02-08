% programme à executer pour resoudre les equation differentiels
clc;clear ALL;close ALL;
% t0=0;
% tf=10;
% delta_t=0.01;
ic=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
% for t=t0:delta_t:tf
% [tt,y]=ode45(@(tt,y)dwESO(y),[t0 t0+delta_t],ic);
% t0=t0+delta_t
% ic=y(end,:)
% end
[T,y]=ode45(@(T,y)Copy_of_dwESO(y),[0 10],ic);
T=0:10/(length(T)-1):10;
figure
subplot(3,1,1)
plot(T,y(:,1),'r') % phi
title('Angle de Roulis')
ylabel('\phi °','fontsize',22)
xlabel('temps s')
legend('sans ESO')
subplot(3,1,2)
plot(T,y(:,4),'b')  %teta
title('Angle de tangage')
ylabel('\theta °','fontsize',22)
xlabel('temps s')
legend('sans ESO')
subplot(3,1,3)
plot(T,y(:,7),'k')  %psi
title('Angle de lacet')
ylabel('\psi °','fontsize',22)
xlabel('temps s')
legend('sans ESO')

figure
xx=subplot(3,1,1);
plot(T,y(:,13),'r')
yy=subplot(3,1,2);
plot(T,y(:,14),'b')
zz=subplot(3,1,3);
plot(T,y(:,15),'k')
title(xx,'L');
title(yy,'M');
title(zz,'N');