function z=myfun4_2(t,x)
%计算m!=0系数u中单独参数
global ringdensity R B ku lamida miu Eh k3 m M N K a b e phif phir;%定义全局变量；
m1=-ringdensity*R;
kmn=-Eh/(1-miu^2)*(2*m*pi/B)^2-ku*R;
c1=2*lamida*sqrt(kmn/m1);

delta2=zeros(1,M);
s2=zeros(1,K);
s3=zeros(1,K);
q2=0;%m!=0路面合激励
for k=1:K
    s3(k)=s3(k)+(phir(k)-phif(k))/K;
    for j=0:M-1
        for i=0:N-1
                if i==0
                    delta2(j+1)=delta2(j+1)+a(j+1,i+1)*(phir(k)^2-phif(k)^2)/2;
                else
                    z1=(cos(i*phir(k))-cos(i*phif(k)))/(-i^2);%中间变量
                    z2=(sin(i*phir(k))-sin(i*phif(k)))/(-i^2);%中间变量
                    delta2(j+1)=delta2(j+1)+a(j+1,i+1)*z1+b(j+1,i+1)*z2;
                end
        end
        if j==0
            z3=cos(2*m*pi*(-1/2+k/K))-cos(2*m*pi*(-1/2+k/K-1/K));
            s2(k)=delta2(j+1)*(-B/2/m/pi)*z3+e(j+1)*(phir(k)-phif(k))*z3;
        elseif j==m
                z6=B/2/K-B/2/(j+m)/pi*sin((j+m)*pi/K)*cos(2*(j+m)*pi*(-1/2+k/K-1/2/K));
                s2(k)=s2(k)+delta2(j+1)*z6+e(j+1)*(phir(k)-phif(k))*z6;
        else
                z7=B/2/(j-m)/pi*sin((j-m)*pi/K)*cos(2*(j-m)*pi*(-1/2+k/K-1/2/K))-B/2/(j+m)/pi*sin((j+m)*pi/K)*cos(2*(j+m)*pi*(-1/2+k/K-1/2/K));
                s2(k)=s2(k)+delta2(j+1)*z7+e(j+1)*(phir(k)-phif(k))*z7;    
        end
    end
    q2=q2+(k3*R/pi/B)*s2(k);%+d3*k3/2/pi*s3(k)
end
%以上为求解路面激励的程序

%下面是求解带束层位移系数的程序        
z(1)=x(2);
z(2)=-c1/m1*x(2)-kmn/m1*x(1)+q2/m1;
z=z';
end