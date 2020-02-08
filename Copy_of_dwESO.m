function ydot = Copy_of_dwESO( y )
kd1=6.4;kd2=6.4;kd3=6.4;
kp1=16;kp2=16;kp3=16;
phi_c=10;teta_c=-10;psi_c=10;
%--------------------------------------------------------------------------
% Ix=12874.8;
% Iy=75673.6;
% Iz=85552.1;
% Ixz=1331.4;
Ix=22682;
Iy=77095;
Iz=95561;
Ixz=1125;
SOM=Ix*Iz-Ixz^2;
c1=((Iy-Iz)*Iz-Ixz^2)/SOM;
c2=((Ix-Iy+Iz)*Ixz)/SOM;
c3=Iz/SOM;
c4=Ixz/SOM;
c5=(Iz-Ix)/Iy;
c6=Ixz/Iy;
c7=1/Iy;
c8=(Ix*(Ix-Iy)+Ixz^2)/SOM;
c9=Ix/SOM;
%--------------------------------------------------------------------------
phi=y(1);
teta=y(4);
p=y(10);
r=y(11);
q=y(12);
%--------------------------------------------------------------------------
B=[c3+c4*cosd(phi)*tand(teta),c7*sind(phi)*tand(teta),c4+c9*cosd(phi)*tand(teta);
    -c4*sind(phi),c7*cosd(phi),-c9*sind(phi);
    c4*cosd(phi)*(1/cosd(teta)),c7*sind(phi)*(1/cosd(teta)),c9*cosd(phi)*(1/cosd(teta))];
phidot=p+tand(teta)*(q*sind(phi)+r*cosd(phi));
tetadot=q*cosd(phi)-r*sind(phi);
A1=((c1*r+c2*p)*q)+(tetadot/cosd(teta)^2)*(q*sind(phi)+r*cosd(phi))+tand(teta)*((c5*p*r-c6*(p^2-r^2))*sind(phi)+q*(phidot*cosd(phi))+((c8*p-c2*r)*q)*cosd(phi)-r*phidot*sind(phi));
A2=(c5*p*r-c6*(p^2-r^2))*cosd(phi)-phidot*sind(phi)*q-(((c8*p-c2*r)*q)*sind(phi)+phidot*cosd(phi)*r);
A3=((c5*p*r-c6*(p^2-r^2))*sind(phi)+phidot*cosd(phi)*q+((c8*p-c2*r)*q)*cosd(phi)-phidot*sind(phi)*r)*cosd(teta)/cosd(teta)^2+(tetadot*sind(teta)*(q*sind(phi)+r*cosd(phi)))/cosd(teta)^2;
A=[A1;A2;A3];
%--------------------------------------------------------------------------
Bdot=[c4*cosd(phi)*(tand(teta)^2 + 1)*tetadot - c4*sind(phi)*tand(teta)*phidot,c7*cosd(phi)*tand(teta)*phidot + c7*sind(phi)*(tand(teta)^2 + 1)*tetadot,c9*cosd(phi)*(tand(teta)^2 + 1)*tetadot - c9*sind(phi)*tand(teta)*phidot;
     -c4*cosd(phi)*phidot,-c7*sind(phi)*phidot,-c9*cosd(phi)*phidot;
     (c4*cosd(phi)*sind(teta)*tetadot)/cosd(teta)^2 - (c4*sind(phi)*phidot)/cosd(teta), (c7*cosd(phi)*phidot)/cosd(teta) + (c7*sind(phi)*sind(teta)*tetadot)/cosd(teta)^2, (c9*cosd(phi)*sind(teta)*tetadot)/cosd(teta)^2 - (c9*sind(phi)*phidot)/cosd(teta)];
%--------------------------------------------------------------------------
pp=inv(B);
L=pp(1,:)*[(kp1*(phi_c-y(1))-kd1*y(2))-A1;((kp2*(teta_c-y(4))-kd2*y(5))-A2);((kp3*(psi_c-y(7))-kd3*y(8))-A3)];
M=pp(2,:)*[(kp1*(phi_c-y(1))-kd1*y(2))-A1;((kp2*(teta_c-y(4))-kd2*y(5))-A2);((kp3*(psi_c-y(7))-kd3*y(8))-A3)];
N=pp(3,:)*[(kp1*(phi_c-y(1))-kd1*y(2))-A1;((kp2*(teta_c-y(4))-kd2*y(5))-A2);((kp3*(psi_c-y(7))-kd3*y(8))-A3)];
pdot=(c1*r+c2*p)*q+c3*L+c4*N;
rdot=(c8*p-c2*r)*q+c4*L+c9*N;
qdot=c5*p*r-c6*(p^2-r^2)+c7*M;

phidotdot=pdot+(tetadot/cosd(teta)^2)*(q*sind(phi)+r*cosd(phi))+tand(teta)*(qdot*sind(phi)+q*(phidot*cosd(phi))+rdot*cosd(phi)-r*phidot*sind(phi));
tetadotdot=qdot*cosd(phi)-phidot*sind(phi)*q-(rdot*sind(phi)+phidot*cosd(phi)*r);

A1dot=(c1*rdot+c2*pdot)*q+qdot*(c1*r+c2*p)+(tetadotdot*cosd(teta)^2+2*tetadot^2*sind(teta)*cosd(teta))/cosd(teta)^4*(q*sind(phi)+r*cosd(phi))+(qdot*sind(phi)+phidot*cosd(phi)*q+rdot*cosd(phi)-phidot*sind(phi)*r)*(tetadot/cosd(teta)^2)+(tetadot/cosd(teta)^2)*((c5*p*r-c6*(p^2-r^2))*sind(phi)+q*(phidot*cosd(phi))+((c8*p-c2*r)*q)*cosd(phi)-r*phidot*sind(phi))+((c5*p*r-c6*(p^2-r^2))*phidot*cosd(phi)+(c5*(pdot*r+rdot*p)-c6*(2*p*pdot-2*r*rdot))*sind(phi)+(qdot*phidot+phidotdot*q)*cosd(phi)-phidot*sind(phi)*q*phidot+(c8*pdot-c2*rdot)*q*cosd(phi)+(qdot*cosd(phi)-phidot*sind(phi)*q)*(c8*p-c2*r)-((rdot*phidot+phidotdot*r)*sind(phi)+phidot*cosd(phi)*r*phidot))*tand(teta);
A2dot=(c5*(pdot*r+rdot*p)-c6*(2*p*pdot-2*r*rdot))*cosd(phi)-phidot*sind(phi)*(c5*p*r-c6*(p^2-r^2))-((phidotdot*q+qdot*phidot)*sind(phi)+phidot*cosd(phi)*phidot*q)-(((c8*pdot-c2*rdot)*q+qdot*(c8*p-c2*r))*sind(phi)+phidot*cosd(phi)*((c8*p-c2*r)*q)+(phidotdot*r+rdot*phidot)*cosd(phi)-phidot*sind(phi)*phidot*r);
A3dot=((((c5*(pdot*r+rdot*p)-c6*(2*p*pdot-2*r*rdot))*sind(phi)+phidot*cosd(phi)*(c5*p*r-c6*(p^2-r^2))+(phidotdot*cosd(phi)-phidot*sind(phi)*phidot)*q+qdot*phidot*cosd(phi)+((c8*pdot-c2*rdot)*q+qdot*(c8*p-c2*r))*cosd(phi)-phidot*sind(phi)*((c8*p-c2*r)*q)-((phidotdot*sind(phi)+phidot*cosd(phi)*phidot)*r+rdot*phidot*sind(phi)))*cosd(teta)-tetadot*sind(teta)*((c5*p*r-c6*(p^2-r^2))*sind(phi)+phidot*cosd(phi)*q+(c5*p*r-c6*(p^2-r^2))*cosd(phi)-phidot*sind(phi)*r))*cosd(teta)^2+2*tetadot*cosd(teta)*sind(teta)*((c5*p*r-c6*(p^2-r^2))*sind(phi)+phidot*cosd(phi)*q+((c8*p-c2*r)*q)*cosd(phi)-phidot*sind(phi)*r)*cosd(teta))/cosd(teta)^4+(((tetadotdot*sind(teta)+tetadot*cosd(teta)*tetadot)*(q*sind(phi)+r*cosd(phi))+(qdot*sind(phi)+phidot*cosd(phi)*q+rdot*cosd(phi)-phidot*sind(phi)*r)*tetadot*sind(teta))*cosd(teta)^2+2*tetadot*cosd(teta)*sind(teta)*(tetadot*sind(teta)*(q*sind(phi)+r*cosd(phi))))/cosd(teta)^4;
Adot=[A1dot;A2dot;A3dot];
%--------------------------------------------------------------------------
u=[L;M;N];
ppdot=[ 0, -(c4*cosd(phi)*phidot)/(c4^2*cosd(phi)^2 + c4^2*sind(phi)^2 - c3*c9*cosd(phi)^2 - c3*c9*sind(phi)^2), (cosd(teta)*(c9*cosd(phi)^2*(tand(teta)^2 + 1)*tetadot - c4*sind(phi)*phidot + c9*sind(phi)^2*(tand(teta)^2 + 1)*tetadot))/(c4^2*cosd(phi)^2 + c4^2*sind(phi)^2 - c3*c9*cosd(phi)^2 - c3*c9*sind(phi)^2) - (sind(teta)*tetadot*(c9*tand(teta)*cosd(phi)^2 + c4*cosd(phi) + c9*tand(teta)*sind(phi)^2))/(c4^2*cosd(phi)^2 + c4^2*sind(phi)^2 - c3*c9*cosd(phi)^2 - c3*c9*sind(phi)^2);
        0,                                                    -(sind(phi)*phidot)/(c7*cosd(phi)^2 + c7*sind(phi)^2),                                                                                                                                                                                                                                                                                             (cosd(phi)*cosd(teta)*phidot)/(c7*cosd(phi)^2 + c7*sind(phi)^2) - (sind(phi)*sind(teta)*tetadot)/(c7*cosd(phi)^2 + c7*sind(phi)^2);
        0,  (c3*cosd(phi)*phidot)/(c4^2*cosd(phi)^2 + c4^2*sind(phi)^2 - c3*c9*cosd(phi)^2 - c3*c9*sind(phi)^2), (sind(teta)*tetadot*(c4*tand(teta)*cosd(phi)^2 + c3*cosd(phi) + c4*tand(teta)*sind(phi)^2))/(c4^2*cosd(phi)^2 + c4^2*sind(phi)^2 - c3*c9*cosd(phi)^2 - c3*c9*sind(phi)^2) - (cosd(teta)*(c4*cosd(phi)^2*(tand(teta)^2 + 1)*tetadot - c3*sind(phi)*phidot + c4*sind(phi)^2*(tand(teta)^2 + 1)*tetadot))/(c4^2*cosd(phi)^2 + c4^2*sind(phi)^2 - c3*c9*cosd(phi)^2 - c3*c9*sind(phi)^2)];
v=[(kp1*(phi_c-y(1))-kd1*y(2));(kp2*(teta_c-y(4))-kd2*y(5));(kp3*(psi_c-y(7))-kd3*y(8))];
vdot=[-kp1*y(2)-kd1*(y(3)+kp1*phi_c);-kp2*y(5)-kd2*(y(6)+kp2*teta_c);-kp3*y(8)-kd3*(y(9)+kp3*psi_c)];
udot=ppdot*(v-A)+pp*(vdot-Adot);
%--------------------------------------------------------------------------
ydot(1)=y(2);   
%ydot(2)=A1+B(1,:)*[L;M;N];
ydot(2)=y(3)+kp1*phi_c;
ydot(3)=A1dot+Bdot(1,:)*u+B(1,:)*udot;
%--------------------------------------------------------------------------
ydot(4)=y(5);   
%ydot(4)=A2+B(2,:)*[L;M;N];  
ydot(5)=y(6)+kp2*teta_c;
ydot(6)=A2dot+Bdot(2,:)*u+B(2,:)*udot;
%--------------------------------------------------------------------------
ydot(7)=y(8);   
%ydot(6)=A3+B(3,:)*[L;M;N];
ydot(8)=y(9)+kp3*psi_c;
ydot(9)=A3dot+Bdot(3,:)*u+B(3,:)*udot;
%--------------------------------------------------------------------------
ydot(10)=(c1*r+c2*p)*q+c3*L+c4*N;  %=pdot
ydot(11)=(c8*p-c2*r)*q+c4*L+c9*N;  %=rdot
ydot(12)=c5*p*r-c6*(p^2-r^2)+c7*M; %=qdot

ydot=[ydot(1);ydot(2);ydot(3);ydot(4);ydot(5);ydot(6);ydot(7);ydot(8);ydot(9);ydot(10);ydot(11);ydot(12);L;M;N];
ydot=ydot(:);
end
