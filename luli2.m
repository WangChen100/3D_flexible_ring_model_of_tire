function [Fx Fy]=luli2(x,qx,N,omega,R,d0)
cx=2.1e6;%胎面纵向总刚度
cy=1070800;%胎面侧向总刚度
cpx=cx/(x(N)-x(1));%胎面纵向单位面积上刚度
cpy=cy/(x(N)-x(1));%胎面侧向单位面积上刚度
vc=omega*(R-d0);
alfa=1.5;
uax=1.05;uay=0.85;n=0.0102;
for s=0:0.02:0.6
    vx=vc*(1+s);
    vsx=20*s;
    vsy=vx*tan(alfa*pi/180);
    usx=uax-n*vsx;%纵侧向摩擦系数
    usy=uay-n*vsy;
    belta=atan(vsy*usy/vsx/usx);
    fx1=sqrt((cpx*s/uax)^2+(cpy*vsy/vc/uay)^2)*(x-x(1));
    fx2=fx1-qx;
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
Fx(t)=cpx*s/2*(xs-x(1))^2+usx*qx2*cos(belta);
Fy(t)=cpy*(vsy/vc)/2*(xs-x(1))^2+usy*qx2*sin(belta);
end
end