function z=myfun1_1(t,x)
%计算m=0,n=0是的系数a、b
global ringdensity omega R B p0 sigmah kw kv ku lamida D k1 k2 k3 M N K n a b phif phir;%定义全局变量；
m1=ringdensity*(n^2+2);
g=2*n^3*ringdensity*omega;
kmn=ringdensity*omega^2*(-n^4+n^2-1)+sigmah/R^2*(n^4+n^2+1)-p0/R+(D/R^4)*(n^6-2*n^4+n^2)+n^2*kw+kv+ku;
c1=2*lamida*sqrt(m1*kmn);

delta1=zeros(1,M);
s1=zeros(1,K);
q1=0;%n=0路面合激励
for k=1:K
    for j=0:M-1
        for i=0:N-1
            if i==0
                delta1(j+1)=delta1(j+1)+(-k2-k3)*a(j+1,i+1)*(phir(k)-phif(k));
            else
                z1=sin(i*phir(k))-sin(i*phif(k));%中间变量
                z2=cos(i*phir(k))-cos(i*phif(k));%中间变量
                delta1(j+1)=delta1(j+1)+(-i^2*k1(k)-k2-k3)/i*(a(j+1,i+1)*z1-b(j+1,i+1)*z2);
            end
        end
        if j==0
            s1(k)=delta1(j+1)*B/K;
        else
            z7=cos((2*j*pi/B)*(-B/2+k*B/K))-cos((2*j*pi/B)*(-B/2+k*B/K-B/K));
            s1(k)=s1(k)+delta1(j+1)*(-B/2/j/pi)*z7;
        end
    end
    q1=q1+(0.5/pi/B)*s1(k);
end
%以上为求解路面激励的程序

%下面是求解带束层位移系数的程序        
z21=x(3);
z22=x(4);
z23=-c1/m1*x(3)-g/m1*x(4)-kmn/m1*x(1)+q1/m1;
z24=-c1/m1*x(4)+g/m1*x(3)-kmn/m1*x(2);
z=[z21;z22;z23;z24];
end