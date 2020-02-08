% programme à executer pour resoudre les equation differentiels
clc;clear ALL;close ALL;
[t,y]=ode45(@(t,y)fslve(y),[0 10],[0 0]);
figure 
plot(y(:,1))
figure 
plot(y(:,2))