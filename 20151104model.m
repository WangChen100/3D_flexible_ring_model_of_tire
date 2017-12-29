clc;
clear all;

global ringdensity omega R B h p0 sigmah sigmahy kw kv ku lamida miu mm I1 I3 tau Eh D G k1 k2 k3 M N K m n a b c d e phif phir Fx Fy Fz T M1 M2 L v1 v2 v3 X;%定义全局变量；
ringdensity=22.4;%rou*h
omega=0;%初始角速度为0
r=0.2032;%轮毂的半径，16英寸
R=0.31;
B=0.165;
h=0.0125;
p0=2.5e5;%胎压
middelR=0.27;sigmahy=p0*(R^2-middelR^2)/2/R/h;%带束层轴向预应力
sigmah=p0*R+ringdensity*R^2*omega^2+2*sigmahy*h*R*tan(pi/6)/B;%带束层周向预应力
kw=5.3e4/B;%胎体径向刚度
kv=7.59e4/B;%胎体切向刚度
ku=7.23e4/B;
lamida=0.2;%阻尼比
miu=0.49;%泊松比
mm=12;%轮毂质量
I1=mm*r^2/2+mm*B^2/12;%有问题！！！！！
I3=mm*r^2;
tau=0.018; %胎面橡胶厚度
E=1142.85;
Eh=14.2857;
D=Eh*h^2/12/(1-miu^2);
G=E/2/(1+miu);


M=5;%模态数15
N=5;%模态数20
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

Fx=4000;d0=0.0135;
Fy=0;
% Fz=2.05*0.3*(omega*(0.328-d0))^2/16/4*9.8;
% T=1.5*Fz*(0.328-d0);
Fz=0;
T=0;
M1=0/4;
M2=0/6;

v1=zeros(1,K);%前角v值
v2=zeros(1,K);%后角v值
v3=zeros(1,K);%后角v值
X=zeros(4,K);%插值参数！


y=-B/2:B/(K-1):B/2;
BB=B*1.3;
k1=(2.54e6/B)*(1-0.3*cos(2*pi*y/BB)+0.3*cos(6*pi*y/BB)/3-0.3*cos(10*pi*y/BB)/5);%胎面径向刚度
k2=6.294e7/B;%胎面切向刚度  有问题！！！
k3=2.656e6/B;%胎面侧向刚度  有问题！！！

phi=-pi/6:2*pi/6/(L-1):pi/6;
x=(R+tau-d0)*tan(phi);

sqacetime=0.2;
sqacetime2=0.01;

tic;
for time=1:1
    time
%     if time==26
%     omega=50;
%     Fz=2.05*0.3*(omega*(0.328-d0))^2/16;
%     T=1.2*Fz*(0.328-d0);
%     end
%     if Fx<1000
%        Fx=500*time;
%    elseif Fx>=1000&&Fx<4000
%        Fx=1000+1000*(time-2);
%     end
% if time>5
%     M1=30;
% end
    for m=0:M-1
        for n=0:N-1
            if m==0
                if n==0
                    [t f11]=ode23(@myfun1_1,0:sqacetime2:sqacetime,[a(m+1,n+1);b(m+1,n+1);va(m+1,n+1);vb(m+1,n+1)]);
                    a2(m+1,n+1)=f11(end,1);
                    b2(m+1,n+1)=f11(end,2);
                    va2(m+1,n+1)=f11(end,3);
                    vb2(m+1,n+1)=f11(end,4);
                elseif n==1
                    [t f12]=ode23(@myfun1_2,0:sqacetime2:sqacetime,[a(m+1,n+1);b(m+1,n+1);ztr1;ztr2;va(m+1,n+1);vb(m+1,n+1);vztr1;vztr2]);
                    a2(m+1,n+1)=f12(end,1);
                    b2(m+1,n+1)=f12(end,2);
                    ztr12=f12(end,3);
                    ztr22=f12(end,4);
                    va2(m+1,n+1)=f12(end,5);
                    vb2(m+1,n+1)=f12(end,6);
                    vztr12=f12(end,7);
                    vztr22=f12(end,8);
                else
                    [t f13]=ode23(@myfun1_3,0:sqacetime2:sqacetime,[a(m+1,n+1);b(m+1,n+1);va(m+1,n+1);vb(m+1,n+1)]);
                    a2(m+1,n+1)=f13(end,1);
                    b2(m+1,n+1)=f13(end,2);
                    va2(m+1,n+1)=f13(end,3);
                    vb2(m+1,n+1)=f13(end,4);
                end
            else
                if n==0
                    [t f21]=ode23(@myfun2_1,0:sqacetime2:sqacetime,[a(m+1,n+1);b(m+1,n+1);va(m+1,n+1);vb(m+1,n+1)]);
                    a2(m+1,n+1)=f21(end,1);
                    b2(m+1,n+1)=f21(end,2);
                    va2(m+1,n+1)=f21(end,3);
                    vb2(m+1,n+1)=f21(end,4);
                else
                    [t f22]=ode23(@myfun2_2,0:sqacetime2:sqacetime,[a(m+1,n+1);b(m+1,n+1);va(m+1,n+1);vb(m+1,n+1)]);
                    a2(m+1,n+1)=f22(end,1);
                    b2(m+1,n+1)=f22(end,2);
                    va2(m+1,n+1)=f22(end,3);
                    vb2(m+1,n+1)=f22(end,4);
                end
            end
            
            if m==0
            if n==0
                [t f31]=ode23(@myfun3_1,0:sqacetime2:sqacetime,[c(1,n+1);d(1,n+1);ztzeta;vc(1,n+1);vd(1,n+1);vztzeta]);
                c2(1,n+1)=f31(end,1);
                d2(1,n+1)=f31(end,2);
                ztzeta2=f31(end,3);
                vc2(1,n+1)=f31(end,4);
                vd2(1,n+1)=f31(end,5);
                vztzeta2=f31(end,6);
            elseif n==1
                [t f32]=ode23(@myfun3_2,0:sqacetime2:sqacetime,[c(1,n+1);d(1,n+1);ztx;ztz;vc(1,n+1);vd(1,n+1);vztx;vztz]);
%                 ff(1+(time-1)*10:11+(time-1)*10,:)=f;
                c2(1,n+1)=f32(end,1);
                d2(1,n+1)=f32(end,2);
                ztx2=f32(end,3);
                ztz2=f32(end,4);
                vc2(1,n+1)=f32(end,5);
                vd2(1,n+1)=f32(end,6);
                vztx2=f32(end,7);
                vztz2=f32(end,8);
            else
                [t f33]=ode23(@myfun3_3,0:sqacetime2:sqacetime,[c(1,n+1);d(1,n+1);vc(1,n+1);vd(1,n+1)]);
                c2(1,n+1)=f33(end,1);
                d2(1,n+1)=f33(end,2);
                vc2(1,n+1)=f33(end,3);
                vd2(1,n+1)=f33(end,4);
            end
            end
        end
        if m==0
            [t f41]=ode23(@myfun4_1,0:sqacetime2:sqacetime,[e(m+1);zty;ve(m+1);vzty]);
            e2(m+1)=f41(end,1);
            zty2=f41(end,2);
            ve2(m+1)=f41(end,3);
            vzty2=f41(end,4);
        else
            [t f42]=ode23(@myfun4_2,0:sqacetime2:sqacetime,[e(m+1);ve(m+1)]);
            e2(m+1)=f42(end,1);
            ve2(m+1)=f42(end,2);
        end
    end
    
    clear f11 f12 f13 f21 f22 f31 f32 f33 f41 f42;
    
   a=a2;%中间变量
   b=b2;
   c=c2;
   d=d2;
   e=e2;
%    va=va2;
%    vb=vb2;
%    vc=vc2;
%    vd=vd2;
%    ve=ve2;
   
   ztx=ztx2;
   ztz=ztz2;
   zty=zty2;
   ztr1=ztr12;
   ztr2=ztr22;
   ztzeta=ztzeta2;
%    vztx=vztx2;
%    vztz=vztz2;
%    vzty=vzty2;
%    vztzeta=vztzeta2;
%    vztr1=vztr12;
%    vztr2=vztr22;

    
   w=zeros(L,K);
   qw=zeros(L,K);
   v=zeros(L,K);
   qv=zeros(L,K);
   u=zeros(L,K);
   qu=zeros(L,K);
   qz=zeros(L,K);
%    Fsx=zeros(101,K);%变滑移率
%    Fsy=zeros(101,K);%变滑移率
%    ms=zeros(101,K);%变滑移率
    Fsx=zeros(L,K);%变滑移率
    Fsy=zeros(L,K);%变滑移率
    ms=zeros(L,K);%变滑移率
   area=0;
   for k=1:K
% %        d1=d0+y(k)*tan(r1);
           for j=0:M-1
               if j==0
                   for i=0:N-1 
                       w(:,k)=(w(:,k)'+(i*a(j+1,i+1)*sin(i*phi)-i*b(j+1,i+1)*cos(i*phi))*y(k))';
                       v(:,k)=(v(:,k)'+(a(j+1,i+1)*cos(i*phi)+b(j+1,i+1)*sin(i*phi))*y(k))';
                       if i==0
                           u(:,k)=(u(:,k)'+a(j+1,i+1)*phi)';
                       else
                           u(:,k)=(u(:,k)'+(a(j+1,i+1)/i*sin(i*phi)-b(j+1,i+1)/i*cos(i*phi)))';
                       end
                   end
                   u(:,k)=u(:,k)+e(j+1)-a(j+1,1)*omega*sqacetime*time;
               else
                   for i=0:N-1 
                       w(:,k)=(w(:,k)'+(i*a(j+1,i+1)*sin(i*phi)-i*b(j+1,i+1)*cos(i*phi))*(-B/2/j/pi)*cos(2*j*pi*y(k)/B))';
                       v(:,k)=(v(:,k)'+(a(j+1,i+1)*cos(i*phi)+b(j+1,i+1)*sin(i*phi))*(-B/2/j/pi)*cos(2*j*pi*y(k)/B))';
                       if i==0
                           u(:,k)=(u(:,k)'+a(j+1,i+1)*phi*sin(2*j*pi*y(k)/B))';
                       else
                           u(:,k)=(u(:,k)'+(a(j+1,i+1)/i*sin(i*phi)-b(j+1,i+1)/i*cos(i*phi))*sin(2*j*pi*y(k)/B))';
                       end
                   end
                   u(:,k)=u(:,k)+(e(j+1)-a(j+1,1)*omega*sqacetime*time)*sin(2*j*pi*y(k)/B);
               end  
           end
           for i=0:N-1
               w(:,k)=(w(:,k)'+i*c(1,i+1)*sin(i*phi)-i*d(1,i+1)*cos(i*phi))';
               v(:,k)=(v(:,k)'+c(1,i+1)*cos(i*phi)+d(1,i+1)*sin(i*phi))';
           end
           qw(:,k)=k1(k)*(w(:,k)'-(R+tau)*(1-cos(phi)))';
           [phif(k) phir(k) z1 z2]=bjtj(phi,qw(:,k));%再编个找出前角和后角的位置 
           qw(qw<0)=0;
           v1(k)=v(z1,k);
           v2(k)=v(z2,k);
           [z5,z6]=min((abs(v(z1:z2,k))));
           z3=z1+z6-1;
%            z3=find(abs(v(:,k)-z5)<1e-8);
           v3(k)=v(z3,k);
           X(:,k)=vnihe(phi(z1),phi(z2),phi(z3),v1(k),v2(k),v3(k));
           z4=X(1,k)*phi(z1:z2).^3+X(2,k)*phi(z1:z2).^2+X(3,k)*phi(z1:z2)+X(4,k);
%            z4=vnihe2(phi,phi(z1),phi(z2),phi(z3),v1(k),v2(k),v3(k),z1,z2,z3);
%            z4=vnihe3(phi(z1:z2),phi(z1),phi(z2),phi(z3),v1(k),v2(k),v3(k));
           qv(z1:z2,k)=k2*(z4'-v(z1:z2,k));
           qv(1:z1,k)=0;
           qv(z2:L,k)=0;
           qu(:,k)=-k3*(-R*u(:,k));
           qz(:,k)=qw(:,k).*cos(phi')+qv(:,k).*sin(phi');
%            [Fsx(:,k),Fsy(:,k),ms(:,k)]=brush1(x(z1:z2),z1,z2,qz(z1:z2,k),omega,R+tau,d0,B,K);
%            [Fsx(:,k),Fsy(:,k),ms(:,k)]=brush2(x(z1:z2),z1,z2,qz(z1:z2,k),omega,R+tau,d0,B,K);%待修改
           [Fsx(z1:z2,k),Fsy(z1:z2,k),ms(z1:z2,k)]=brush3(x(z1:z2),z1,z2,qz(z1:z2,k),omega,R+tau,d0,B,K);
           area=area+(x(z2)-x(z1))*B/K;    
   end
end
toc
% ms
area%接地面积
qzsum=sum(qz,2);
sum(qzsum(1:(L+1)/2))/sum(qzsum((L+1)/2:L))%纵向非对称指数
sum(qz(:,8).*x')/sum(qz(:,8))%中面滚动阻力距
meanms=sum(ms);

% save 20151125
% figure(3)
% meshz(y,x,qv);
% colorbar('vertical')
% axis tight    
% xlabel('轮胎宽度方向','fontsize',8)
% ylabel('轮胎圆周方向','fontsize',8)
% zlabel('路面径向激励分布N/m2','fontsize',8)
% hold on;
figure(4)
meshz(y,x,qz);
colorbar('vertical')
axis([-0.0825 0.0825 -0.12 0.12]);  
xlabel('轮胎宽度方向/m','FontName','宋体','fontsize',14)
ylabel('轮胎圆周方向/m','FontName','宋体','fontsize',14)
zlabel('轮胎接地压力分布N/m^2','FontName','宋体','fontsize',14)
s=0:0.01:1;
% figure(5)
% meshz(y,s,ms);
% colorbar('vertical')
% xlabel('轮胎宽度方向','FontName','宋体','fontsize',14)
% ylabel('滑移率','FontName','宋体','fontsize',14)
% zlabel('胎面磨损','FontName','宋体','fontsize',14)

figure(6)
surf(y+0.0057,x,qz);
colorbar('vertical')
axis([-0.09 0.09 -0.1 0.1]);  
xlabel('轮胎宽度方向/m','FontName','宋体','fontsize',14)
ylabel('轮胎圆周方向/m','FontName','宋体','fontsize',14)
zlabel('轮胎接地压力分布N/m^2','FontName','宋体','fontsize',14)
figure(7)
plot(y,maxms,'-k','linewidth',2);
axis tight
xlabel('轮胎宽度方向/m','FontName','宋体','fontsize',14)
ylabel('轮胎磨损/Kg/m','FontName','宋体','fontsize',14)
figure(8)
surf(y,x,ms);
colorbar('vertical')
axis([-0.0825 0.0825 -0.12 0.12]);
set(gcf,'color','white')
grid off
xlabel('轮胎宽度方向/m','FontName','宋体','fontsize',14)
ylabel('轮胎圆周方向/m','FontName','宋体','fontsize',14)
zlabel('胎面磨损量/Kg/m^2','FontName','宋体','fontsize',14)

figure(1)
plot(x(9:33),R*u(9:33,1),'--*k','linewidth',2);
axis tight
xlabel('轮胎纵向/m','FontName','宋体','fontsize',14)
ylabel('接地印迹侧向偏移量/m','FontName','宋体','fontsize',14)
hold on

figure(1)
plot3(y(1)*ones(size(x(7:35))),x(7:35),2.2*qz(7:35,1),'-k','linewidth',3);
hold on
plot3(y(4)*ones(size(x(7:35))),x(7:35),2.2*qz(7:35,4),'-k','linewidth',3);
hold on
plot3(y(8)*ones(size(x(7:35))),x(7:35),2.2*qz(7:35,8),'-k','linewidth',3);
hold on
plot3(y(12)*ones(size(x(7:35))),x(7:35),2.2*qz(7:35,12),'-k','linewidth',3);
hold on
plot3(y(15)*ones(size(x(7:35))),x(7:35),2.2*qz(7:35,15),'-k','linewidth',3);
hold on
grid on
colorbar('vertical')
axis tight    
xlabel('轮胎宽度方向','fontsize',14)
ylabel('轮胎圆周方向','fontsize',14)
zlabel('路面径向激励分布(pa)','fontsize',14)


figure(2)
plot(y,meanms,'-^k','linewidth',3);
hold on
axis tight
xlabel('轮胎侧向（m）','FontName','宋体','fontsize',14)
ylabel('单位面积轮胎磨损量（mg/mm^2）','FontName','宋体','fontsize',14)

% figure(1)
% meshz(y,phi,qw);
% colorbar('vertical')
% axis tight    
% xlabel('轮胎宽度方向','fontsize',8)
% ylabel('轮胎圆周方向','fontsize',8)
% zlabel('路面径向激励分布N/m2','fontsize',8)
% figure(2)
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
