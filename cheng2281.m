clc;
clear all; 
omega=100;
lamda=0.15;         
R=0.3;
b=0.14;
tau=0.005;                 
ringdensity=3.15; 
kw=6.3*10^5;
kv=1.89*10^5;
p0=1.2*10^5; 
EI=2;

sigmaA=p0*b*R+ringdensity*R^2*omega^2; 
kk2=4.77e7;           %second spring stifness
 
N=40;
j=1:N;
kn=(EI*j.^2/R^4+sigmaA/R^2).*(1-j.^2).^2-p0*b*(1-j.^2)/R+kv+kw*j.^2-ringdensity*(1+j.^2).*omega^2;
mn=ringdensity*(1+j.^2);
gn=-4*ringdensity.*j*omega;
cn=2*lamda*sqrt(mn.*kn); 
rn=(1./j).*atan((cn.*j*omega)./(kn-mn.*(j*omega).^2-gn.*j*omega));
An=-(j.^2/pi)./sqrt((kn-mn.*(j*omega).^2-gn.*j*omega).^2+(cn.*j*omega).^2);
 
 
Phif=-15.1488/180*pi;
Phir=12.4795/180*pi;

i=1;
while i<=N;
    j=1;
    while j<=N;
        z2=sin(j*(Phir-rn(j)));%z2-z7均为中间变量
        z3=cos(j*(Phir-rn(j)));
        z4=(i-j)*Phir+j*rn(j);
        z5=(i-j)*Phif+j*rn(j);
        z6=(i+j)*Phir-j*rn(j);
        z7=(i+j)*Phif-j*rn(j);
        if i~=j
            AA(i,j)=An(j)*((-z3)*(sin(i*Phir)-sin(i*Phif))/i+(sin(z4)-sin(z5))/2/(i-j)+(sin(z6)-sin(z7))/2/(i+j));
            BB(i,j)=An(j)*((-z2)*(sin(i*Phir)-sin(i*Phif))/i+(cos(z4)-cos(z5))/2/(i-j)-(cos(z6)-cos(z7))/2/(i+j));
            CC(i,j)=An(j)*(z3*(cos(i*Phir)-cos(i*Phif))/i-(cos(z4)-cos(z5))/2/(i-j)-(cos(z6)-cos(z7))/2/(i+j));
            DD(i,j)=An(j)*(z2*(cos(i*Phir)-cos(i*Phif))/i+(sin(z4)-sin(z5))/2/(i-j)-(sin(z6)-sin(z7))/2/(i+j));
        else
            AA(i,j)=An(j)*((-z3)*(sin(i*Phir)-sin(i*Phif))/i+cos(j*rn(j))*(Phir-Phif)/2+(sin(z6)-sin(z7))/2/(i+j))-1/kk2;
            BB(i,j)=An(j)*((-z2)*(sin(i*Phir)-sin(i*Phif))/i-sin(j*rn(j))*(Phir-Phif)/2-(cos(z6)-cos(z7))/2/(i+j));
            CC(i,j)=An(j)*(z3*(cos(i*Phir)-cos(i*Phif))/i+sin(j*rn(j))*(Phir-Phif)/2-(cos(z6)-cos(z7))/2/(i+j));
            DD(i,j)=An(j)*(z2*(cos(i*Phir)-cos(i*Phif))/i+cos(j*rn(j))*(Phir-Phif)/2-(sin(z6)-sin(z7))/2/(i+j))-1/kk2;
        end
        j=j+1;
    end
    z8=(sin(i*Phir)-sin(i*Phif))/i;
    z9=(cos(i*Phir)-cos(i*Phif))/i;
    z10=(sin((i-1)*Phir)-sin((i-1)*Phif))/2/(i-1);
    z11=(sin((i+1)*Phir)-sin((i+1)*Phif))/2/(i+1);
    z12=(cos((i-1)*Phir)-cos((i-1)*Phif))/2/(i-1);
    z13=(cos((i+1)*Phir)-cos((i+1)*Phif))/2/(i+1);
    if i~=1
        EE(i,1)=-(R+tau)*(-cos(Phir)*z8+z10+z11);
        FF(i,1)=(R+tau)*(-cos(Phir)*z9+z12+z13);
    else
        EE(i,1)=-(R+tau)*(-cos(Phir)*z8+(Phir-Phif)/2+z11);
        FF(i,1)=(R+tau)*(-cos(Phir)*z9+z13);
    end
    i=i+1;
end  %至此矩阵ABCDEF构建完成
 
GG=[AA BB;CC DD];   %式（5.12）的大矩阵
HH=[EE;FF];
XX=inv(GG)*HH;
 
n=50 ;
Phi=Phif:(Phir-Phif)/(n-1):Phir;%求胎面径向位移的求解
sum1=zeros(1,n);
 for i=1:n
     for j=1:N
        z14=cos(j*(Phi(i)-rn(j)));
        z15=sin(j*(Phi(i)-rn(j)));
        sum1(i)=sum1(i)+(An(j)*(XX(j)*z14+XX(N+j)*z15));
     end
 end
fx1=(R+tau)*(1-cos(Phi(1)))-sum1(1);
fx2=kk2*(sum1-(R+tau)*(1-cos(Phi)))+kk2*fx1;
fx3=fx2.*cos(Phi(j))
R*XX(1)
%R*XX(N+1)
%fx4(i)=fx2(i)*sin(Phi(i))

figure(1)
plot(sin(Phi)*(R+tau),fx2,'k-','linewidth',2.5)
hold on;
xlabel('接地长度(m)','fontsize',13)
ylabel('接地压力（N/m）','fontsize',13)

[Fx,Fy]=luli(sin(Phi)*(R+tau),fx2,n,omega,R,fx1);
[num1,txt,rawdata]=xlsread('C:\Users\fan\Desktop\WangC\M文件\newlongitudinal.xls');
[num2,txt,rawdata]=xlsread('C:\Users\fan\Desktop\WangC\M文件\newlateral.xls');
s=0:0.02:0.6;

figure(2)
plot(s,Fx,'LineWidth',3)
hold on
plot(num1(:,1),num1(:,2),'*')
hold on
legend('本文计算结果','实验结果');
xlabel('滑移率','fontsize',14);
ylabel('纵向力(N)','fontsize',14);
figure(3)
plot(s,Fy,'LineWidth',3)
hold on
plot(num2(:,1),num2(:,2),'*')
legend('本文计算结果','实验结果');
xlabel('滑移率','fontsize',14);
ylabel('侧向力(N)','fontsize',14);


%Phif*(R+tau)
%Phir*(R+tau)
(-Phif*(R+tau)+Phir*(R+tau))/2
grid off
x=sin(Phi)*(R+tau);qx=fx2;N=n;d0=fx1;