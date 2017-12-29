function [Fx Fy]=luli(x,qx,N,omega,R,d0)
cx=2.1e6;%胎面纵向总刚度
cy=1010800;%胎面侧向总刚度
cpx=cx/(x(N)-x(1));%胎面纵向单位面积上刚度
cpy=cy/(x(N)-x(1));%胎面侧向单位面积上刚度
vc=omega*(R-d0);
alfa=4;
tanalfa=tan(alfa*pi/180); %弧度
uax=1;uay=0.89;n=0.0035;
for s=0:0.02:0.6
    aa=cpx*s/(1-s)/uax;
    bb=cpy*tanalfa/(1-s)/uay;
    fx1=sqrt(aa^2+bb^2)*(x-x(1));
    fx2=fx1-qx;
    vx=vc/(1-s);
    vsx=vx*s;
    vsy=vx*tanalfa;
    usx=uax-n*vsx;%纵侧向摩擦系数
    usy=uay-n*vsy;
    if fx2(2)>0
        xs=x(1);
        qx2=sum(qx(1:N))*(x(2)-x(1));
    else
        for i=2:N-1
            if fx2(i)*fx2(i+1)<=0
                xs=(x(i)+x(i+1))/2;
                qx2=sum(qx(i+1:N))*(x(2)-x(1));
                break
            end
        end
    end
t=round(s/0.02+1);
cc=usx*s;
dd=usy*tanalfa;
Fx(t)=cpx*s/2/(1-s)*(xs-x(1))^2+usx*qx2*cc./sqrt(cc.^2+dd.^2);
Fy(t)=cpy*tanalfa/2/(1-s)*(xs-x(1))^2+usy*qx2*dd./sqrt(cc.^2+dd.^2);
end
end