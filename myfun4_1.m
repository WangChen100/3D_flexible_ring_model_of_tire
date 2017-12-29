function z=myfun4_1(t,x)
%计算m=0系数u中单独参数、y
global ringdensity R B ku lamida miu mm Eh k3 m M N K a b e phif phir Fy;%定义全局变量；
m1=-ringdensity*R;
kmn=(-Eh/(1-miu^2)*(2*m*pi/B)^2-ku*R)*2;
c1=2*lamida*sqrt(kmn/m1);
ma=-mm/(2*pi*B*R^2);
ka=-ku;
kc=-ku/R;
ca=2*lamida*sqrt(kc/ma);


delta1=zeros(1,M);
s1=zeros(1,K);
s3=zeros(1,K);
q1=0;%m=0路面合激励
for k=1:K
    s3(k)=s3(k)+(phir(k)-phif(k))/K;
    for j=0:M-1
        for i=0:N-1
            if i==0
                    delta1(j+1)=delta1(j+1)+a(j+1,i+1)*(phir(k)^2-phif(k)^2)/2;   
            else
                    z1=(cos(i*phir(k))-cos(i*phif(k)))/(-i^2);%中间变量
                    z2=(sin(i*phir(k))-sin(i*phif(k)))/(-i^2);%中间变量
                    delta1(j+1)=delta1(j+1)+a(j+1,i+1)*z1+b(j+1,i+1)*z2;    
            end
        end
        if j==0
            s1(k)=delta1(j+1)*(B/K)+e(j+1)*(phir(k)-phif(k))*B/K;
        else
            z5=cos(2*j*pi*(-1/2+k/K))-cos(2*j*pi*(-1/2+k/K-1/K));
            s1(k)=s1(k)+delta1(j+1)*(-B/2/j/pi)*z5+e(j+1)*(phir(k)-phif(k))*(-B/2/j/pi)*z5;
        end
    end
    q1=q1+(k3*R/2/pi/B)*s1(k);%+d3*k3/2/pi*s3(k)
end
%以上为求解路面激励的程序

%下面是求解带束层位移系数的程序        
z11=x(3);
z12=x(4);
z13=-c1*x(3)-kmn/m1*x(1)-ka/m1*x(2)+q1/m1;
z14=-ca*x(4)-ka/ma*x(1)-kc/ma*x(2)-Fy/(2*pi*B*R^2)/ma;
z=[z11;z12;z13;z14];
end