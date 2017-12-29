function [Fsx,Fsy,ms]=brush3(x,z1,z2,qz,omega,R,d0,B,K)%纵侧向三维磨损
Fsx=zeros(z2-z1+1,1);
Fsy=zeros(z2-z1+1,1);
if z2~=z1
    cx=2.1e6;%胎面纵向总刚度
    cy=1.0108e5;%胎面侧向总刚度
    cpx=cx/(x(end)-x(1)+eps);%胎面纵向单位面积上刚度
    cpy=cy/(x(end)-x(1)+eps);%胎面侧向单位面积上刚度
    vc=omega*(R-d0);
    alfa=4;
    tanalfa=tan(alfa*pi/180); %弧度
    uax=1;uay=0.89;n=0.0035;
    s=0.03;%滑移率
    aa=cpx*s/uax;
    bb=cpy*tanalfa*(1+s)/uay;
    fx1=sqrt(aa^2+bb^2)*(x-x(1));
    fx2=fx1/K-qz'*B/K;
    vx=vc*(1+s);
    vsx=vc*s;
    vsy=vx*tanalfa;
    usx=uax-n*vsx;%纵侧向摩擦系数
    usy=uay-n*vsy;
    N=length(x);
    ii=0;
    if fx2(2)>0
        xs=x(1);
        ii=1;
%         qx2=sum(qz)*(x(2)-x(1));
    else
        for i=2:N-1
            if fx2(i)*fx2(i+1)<=0
                xs=(x(i)+x(i+1))/2;
                ii=i;
%               qx2=sum(qz(i+1:N))*(x(2)-x(1));
                break
            end
        end
    end
cc=usx*s;
dd=usy*(1+s)*tanalfa;
Fsx(ii+1:N)=usx*qz(ii+1:N)*B/K*(x(2)-x(1))*cc/sqrt(cc^2+dd^2);
Fsy(ii+1:N)=usy*qz(ii+1:N)*B/K*(x(2)-x(1))*dd/sqrt(cc^2+dd^2);
ww=sqrt((Fsx*vsx).^2+(Fsy*vsy).^2);
c1=0.0001;c2=1.4339;
ms=c1*((ww/B*K/(x(2)-x(1))/1e6).^c2);
[i1,j1,value1]=find(qz);
[i2,j2,value2]=find(ms==0);
i3=intersect(i1,i2);
j3=intersect(j1,j2);
ms(i3,j3)=1e-9;
% Fsx=usx*qx2*cc/sqrt(cc^2+dd^2);
% Fsy=usy*qx2*dd/sqrt(cc^2+dd^2);
end
end