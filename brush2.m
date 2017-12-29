function [Fsx,Fsy]=brush2(x,z1,z2,qz,omega,R,d0)%定滑移率单点最大磨损
if z2==z1
    Fsx=0;
    Fsy=0;
else
    cx=2.1e6;%胎面纵向总刚度
    cy=1010800;%胎面侧向总刚度
    cpx=cx/(x(end)-x(1)+eps);%胎面纵向单位面积上刚度
    cpy=cy/(x(end)-x(1)+eps);%胎面侧向单位面积上刚度
    vc=omega*(R-d0);
    alfa=4;
    tanalfa=tan(alfa*pi/180); %弧度
    uax=1;uay=0.89;n=0.0035;
    s=0.1;%滑移率
    aa=cpx*s/uax;
    bb=cpy*tanalfa*(1+s)/uay;
    fx1=sqrt(aa^2+bb^2)*(x-x(1));
    fx2=fx1-qz;
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
Fsx=cpx*s/2*(xs-x(1))^2+usx*qx2*cc/sqrt(cc^2+dd^2);
Fsy=cpy*tanalfa*(1+s)/2*(xs-x(1))^2+usy*qx2*dd/sqrt(cc^2+dd^2);
end
end