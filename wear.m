function [ms]=wear(Fsx,Fsy,B,K,omega,R,d0)
vc=omega*(R-d0);
s=0.03;%������
alfa=4;
tanalfa=tan(alfa*pi/180); %����
vsx=vc*s;
vx=vc*(1+s);
vsy=vx*tanalfa;

w=B/(K-1)*sqrt((sum(Fx)*vsx)^2+(sum(Fy)*vsy)^2);
c1=0.0001;c2=1.4339;
ms=c1*(w^c2);
end