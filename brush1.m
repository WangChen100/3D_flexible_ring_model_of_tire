function [Fsx,Fsy,ms]=brush1(x,z1,z2,qz,omega,R,d0,B,K)%随滑移率变化的磨损
   Fsx=zeros(101,1);
   Fsy=zeros(101,1);
   ms=zeros(101,1);
if z2~=z1
    for j=0:100%滑移率
        s=0.01*j;
    cx=2.1e6;%胎面纵向总刚度
    cy=1010800;%胎面侧向总刚度
    cpx=cx/(x(end)-x(1)+eps);%胎面纵向单位长度上刚度
    cpy=cy/(x(end)-x(1)+eps);%胎面侧向单位长度上刚度
    vc=omega*(R-d0);
    alfa=4;
    tanalfa=tan(alfa*pi/180); %弧度
    uax=1;uay=0.89;n=0.0035;
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
    if fx2(2)>0
        xs=x(1);
        qx2=sum(qz)*(x(2)-x(1));
    else
        for i=2:N-1
            if fx2(i)*fx2(i+1)<=0
                xs=(x(i)+x(i+1))/2;
                qx2=sum(qz(i+1:N))*(x(2)-x(1));
                break
            end
        end
    end
cc=usx*s;
dd=usy*(1+s)*tanalfa;
fx=usx*qx2*B/K*cc/sqrt(cc^2+dd^2);%cpx*s/2*(xs-x(1))^2+
fy=usy*qx2*B/K*dd/sqrt(cc^2+dd^2);%cpy*tanalfa*(1+s)/2*(xs-x(1))^2+
Fsx(j+1,1)=fx;
Fsy(j+1,1)=fy;

ww=sqrt((Fsx(j+1,1)*vsx)^2+(Fsy(j+1,1)*vsy)^2);
c1=0.0001;c2=1.4339;
ms(j+1,1)=c1*(ww.^c2);
    end
   end
end