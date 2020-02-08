% programme à executer pour resoudre les equation differentiels
clc;clear ALL;close ALL;

[t,x]=ode45(@(t,x)controller2(x),[0 10],[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]);

t=0:10/(length(t)-1):10;
figure
subplot(3,1,1)
plot(t,x(:,1),'r') % phi
title('Angle de Roulis')
ylabel('\phi °','fontsize',22)
xlabel('temps s')
legend('avec ESO')
subplot(3,1,2)
plot(t,x(:,4),'b')  %teta
title('Angle de tangage')
ylabel('\theta °','fontsize',22)
xlabel('temps s')
legend('avec ESO')
subplot(3,1,3)
plot(t,x(:,7),'k')  %psi
title('Angle de lacet')
ylabel('\psi °','fontsize',22)
xlabel('temps s')
legend('avec ESO')

figure
aa=subplot(3,1,1);
plot(t,x(:,2),'r')
bb=subplot(3,1,2);
plot(t,x(:,5),'b')
cc=subplot(3,1,3);
plot(t,x(:,8),'k')
title(aa,'\phi^.','fontsize',22);
title(bb,'\theta^.','fontsize',22);
title(cc,'\psi^.','fontsize',22);

figure
a=subplot(3,1,1);
plot(t,x(:,19),'r')
b=subplot(3,1,2);
plot(t,x(:,20),'b')
c=subplot(3,1,3);
plot(t,x(:,21),'k')
title(a,'p','fontsize',22);
title(b,'r','fontsize',22);
title(c,'q','fontsize',22);

kd1=6.4;kd2=6.4;kd3=6.4;
kp1=16;kp2=16;kp3=16;
Delta_u1=x(:,12)+kp1*x(:,1)+kd1*x(:,2);
Delta_u2=x(:,15)+kp2*x(:,4)+kd2*x(:,5);
Delta_u3=x(:,18)+kp3*x(:,7)+kd3*x(:,8);
figure
xx=subplot(3,1,1);
plot(Delta_u1)
yy=subplot(3,1,2);
plot(Delta_u2)
zz=subplot(3,1,3);
plot(Delta_u3)
title(xx,'Delta_u_{1}');
title(yy,'Delta_u_{2}');
title(zz,'Delta_u_{3}');