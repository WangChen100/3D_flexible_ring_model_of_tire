clc;
clear all;

global ringdensity omega R B h p0 sigmah sigmahy kw kv ku lamida miu mm I1 I3 tau Eh D G k1 k2 k3 M N K m n a b c d e phif phir Fx Fy Fz T M1 M2 L v1 v2 v3 X;%定义全局变量；
ringdensity=22.5;%rou*h
omega=0;%初始角速度为0
r=0.2032;%轮毂的半径，16英寸
R=0.3;
B=0.14;
h=0.0125;
p0=1.2e5;%胎压
middelR=0.27;sigmahy=p0*(R^2-middelR^2)/2/R/h;%带束层轴向预应力
sigmah=p0*R+ringdensity*R^2*omega^2+2*sigmahy*h*R*tan(pi/6)/B;%带束层周向预应力
kw=3.3e5/B;%胎体径向刚度
kv=1.89e5/B;%胎体切向刚度
ku=1.23e5/B;
lamida=0.05;%阻尼比
miu=0.49;%泊松比
mm=12;%轮毂质量
I1=mm*B^2/12;%有问题！！！！！
I3=mm*r^2;
tau=0.005; %胎面橡胶厚度
E=1142.85;
Eh=14.2857;
D=Eh*h^2/12/(1-miu^2);
G=E/2/(1+miu);

k1=1.794e6/B;%胎面径向刚度
k2=8.794e7/B;%胎面切向刚度  有问题！！！
k3=2.656e6/B;%胎面侧向刚度  有问题！！！


M=15;%模态数
N=20;%模态数
K=15;%侧向切片数
L=41;%纵向切片数为L-1=40

a=zeros(M,N);%带束层单元位移系数的初始值
b=zeros(M,N);
c=zeros(1,N);
d=zeros(1,N);
e=zeros(M,1);
va=zeros(M,N);%带束层单元位移系数的初始速度值
vb=zeros(M,N);
vc=zeros(M,N);
vd=zeros(M,N);
ve=zeros(M,N);
ztx=0;
zty=0;
ztz=0;
ztzeta=0;
ztr1=0;
ztr2=0;
vztx=0;%轴头自由度的初始速度值
vzty=0;
vztz=0;
vztzeta=0;
vztr1=0;
vztr2=0;

a2=zeros(M,N);%中间变量
b2=zeros(M,N);
c2=zeros(1,N);
d2=zeros(1,N);
e2=zeros(M,1);
va2=zeros(M,N);
vb2=zeros(M,N);
vc2=zeros(M,N);
vd2=zeros(M,N);
ve2=zeros(M,N);
ztx2=0;
zty2=0;
ztz2=0;
ztzeta2=0;
ztr12=0;
ztr22=0;
vztx2=0;
vzty2=0;
vztz2=0;
vztzeta2=0;
vztr12=0;
vztr22=0;

phif=zeros(1,K);
phir=zeros(1,K);

Fx=0;
Fy=0;
Fz=0;
T=0;
M1=0;
M2=0;

v1=zeros(1,K);%前角v值
v2=zeros(1,K);%后角v值
v3=zeros(1,K);%后角v值
X=zeros(4,K);%插值参数！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！

y=-B/2:B/K:B/2;
sqacetime=0.01;
sqacetime2=0.001;

for time=1:30
    time
   if Fx<1000
       Fx=500*time;
   elseif Fx>=1000&&Fx<4000
       Fx=1000+1000*(time-2);
   end       
    for m=0:M-1
        for n=0:N-1            
            if m==0
            if n==0
                [t f]=ode23(@myfun3_1,0:sqacetime2:sqacetime,[c(1,n+1);d(1,n+1);ztzeta;vc(1,n+1);vd(1,n+1);vztzeta]);
                c2(1,n+1)=f(end,1);
                d2(1,n+1)=f(end,2);
                ztzeta2=f(end,3);
                vc2(1,n+1)=f(end,4);
                vd2(1,n+1)=f(end,5);
                vztzeta2=f(end,6);
            elseif n==1
                [t f]=ode23(@myfun3_2,0:sqacetime2:sqacetime,[c(1,n+1);d(1,n+1);ztx;ztz;vc(1,n+1);vd(1,n+1);vztx;vztz]);
                ff(1+(time-1)*10:11+(time-1)*10,:)=f;
                c2(1,n+1)=f(end,1);
                d2(1,n+1)=f(end,2);
                ztx2=f(end,3);
                ztz2=f(end,4);
                vc2(1,n+1)=f(end,5);
                vd2(1,n+1)=f(end,6);
                vztx2=f(end,7);
                vztz2=f(end,8);
            else
                [t f]=ode23(@myfun3_3,0:sqacetime2:sqacetime,[c(1,n+1);d(1,n+1);vc(1,n+1);vd(1,n+1)]);
                c2(1,n+1)=f(end,1);
                d2(1,n+1)=f(end,2);
                vc2(1,n+1)=f(end,3);
                vd2(1,n+1)=f(end,4);
            end
            end
        end
    end
    
   c=c2;
   d=d2;
   ztx=ztx2;
   ztz=ztz2;
   ztzeta=ztzeta2;

   vc=vc2;
   vd=vd2;
   vztx=vztx2;
   vztz=vztz2;
   vzty=vzty2;
   vztzeta=vztzeta2;


    phi=-pi/9:2*pi/9/(L-1):pi/9;
    %omega=vztzeta;%待修改
    w=zeros(L,K);
    qw=zeros(L,K);
    v=zeros(L,K);
    qv=zeros(L,K);
    u=zeros(L,K);
    qu=zeros(L,K);
 
   for k=1:K
%   d1=d0+y(k)*tan(r1);
           for j=0:M-1
               if j==0
                   for i=0:N-1 
                       w(:,k)=(w(:,k)'+(i*a(j+1,i+1)*sin(i*phi)-i*b(j+1,i+1)*cos(i*phi))*y(k))';
                       v(:,k)=(v(:,k)'+(a(j+1,i+1)*cos(i*phi)+b(j+1,i+1)*sin(i*phi))*y(k))';
                       if i==0
                           u(:,k)=a(j+1,i+1)*phi;
                       else
                           u(:,k)=(u(:,k)'+(a(j+1,i+1)/i*sin(i*phi)-b(j+1,i+1)/i*cos(i*phi)))';
                       end
                   end
                   u(:,k)=u(:,k)+e(j+1);
               else
                   for i=0:N-1 
                       w(:,k)=(w(:,k)'+(i*a(j+1,i+1)*sin(i*phi)-i*b(j+1,i+1)*cos(i*phi))*(-B/2/j/pi)*cos(2*j*pi*y(k)/B))';
                       v(:,k)=(v(:,k)'+(a(j+1,i+1)*cos(i*phi)+b(j+1,i+1)*sin(i*phi))*(-B/2/j/pi)*cos(2*j*pi*y(k)/B))';
                       if i==0
                           u(:,k)=a(j+1,i+1)*phi*sin(2*j*pi*y(k)/B);
                       else
                           u(:,k)=(u(:,k)'+(a(j+1,i+1)/i*sin(i*phi)-b(j+1,i+1)/i*cos(i*phi))*sin(2*j*pi*y(k)/B))';
                       end
                   end
                   u(:,k)=u(:,k)+e(j+1)*sin(2*j*pi*y(k)/B);
               end  
           end
           for i=0:N-1
               w(:,k)=(w(:,k)'+i*c(1,i+1)*sin(i*phi)-i*d(1,i+1)*cos(i*phi))';
               v(:,k)=(v(:,k)'+c(1,i+1)*cos(i*phi)+d(1,i+1)*sin(i*phi))';
           end
           qw(:,k)=k1*(w(:,k)'-(R+tau)*(1-cos(phi)))';
           [phif(k) phir(k) z1 z2]=bjtj(phi,qw(:,k));%再编个找出前角和后角的位置 
           qw(qw<0)=0;
           v1(k)=v(z1,k);
           v2(k)=v(z2,k);
           z3=find(abs(v(:,k))<1e-10);
           v3(k)=v(z3,k);
           X(:,k)=vnihe(phi(z1),phi(z2),phi(z3),v1(k),v2(k),v3(k));
           z4=X(1,k)*phi(z1:z2).^3+X(2,k)*phi(z1:z2).^2+X(3,k)*phi(z1:z2)+X(4,k);
           qv(z1:z2,k)=k2*(z4'-v(z1:z2,k));
           qv(1:z1,k)=0;
           qv(z2:41,k)=0;
           qu(:,k)=-k3*(-R*u(:,k));
           qz(:,k)=qw(:,k).*cos(phi')+qv(:,k).*sin(phi');
   end
   
%    figure(1)
% meshz(y,phi,qw);
% colorbar('vertical')
% axis tight    
% xlabel('轮胎宽度方向','fontsize',8)
% ylabel('轮胎圆周方向','fontsize',8)
% zlabel('路面径向激励分布N/m2','fontsize',8)
%    figure(2)
% meshz(y,phi,qv);
% colorbar('vertical')
% axis tight    
% xlabel('轮胎宽度方向','fontsize',8)
% ylabel('轮胎圆周方向','fontsize',8)
% zlabel('路面切向激励分布N/m2','fontsize',8)
% figure(3)
% meshz(y,phi,qu);
% colorbar('vertical')
% axis tight    
% xlabel('轮胎宽度方向','fontsize',8)
% ylabel('轮胎圆周方向','fontsize',8)
% zlabel('路面侧向激励分布N/m2','fontsize',8)
% figure(4)
% meshz(y,phi,qz);
% colorbar('vertical')
% axis tight    
% xlabel('轮胎宽度方向','fontsize',8)
% ylabel('轮胎圆周方向','fontsize',8)
% zlabel('路面垂向压力分布N/m2','fontsize',8)

end