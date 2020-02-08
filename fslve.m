function ydot = fslve(y)
kd1=6.4;
kp1=16;
phi_c=10;
ydot(1)=y(2);
ydot(2)=kp1*(phi_c-y(1))-kd1*y(2);
ydot=[ydot(1);ydot(2)];
ydot=ydot(:);
end
