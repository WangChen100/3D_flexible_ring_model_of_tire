function z=myfun3_1(t,x)
%计算n=0系数c、d
global ringdensity omega R B p0 sigmah kw kv lamida tau D I3 k1 k2 n M N K a b c d phif phir T v1 v2 X;%定义全局变量；
m1=ringdensity*(n^2+1);
g=2*(n^3-n)*ringdensity*omega;
kmn=-ringdensity*omega^2*(n^2-1)^2+p0/R*(n^2-1)+(D/R^4)*(n^6-2*n^4+n^2)+(sigmah/R^2)*(n^4+1)+n^2*kw+kv;
ma=I3/(2*pi*B*R);
ka=-kv*R;
kb=kv*R^2;
c1=2*lamida*sqrt(m1*kmn);
ca=2*lamida*sqrt(ma*kb);

delta11=zeros(1,M);
delta12=zeros(1,K);
s1=zeros(1,K);
q1=0;%n=0路面合激励
for k=1:K
    for j=0:M-1
        for i=0:N-1
            if i==0
                    delta11(j+1)=delta11(j+1)-k2*a(j+1,i+1)*(phir(k)-phif(k));
            else
                    z1=(sin(i*phir(k))-sin(i*phif(k)))/i;%中间变量
                    z2=(cos(i*phir(k))-cos(i*phif(k)))/i;%中间变量
                    delta11(j+1)=delta11(j+1)+(-i^2*k1(k)-k2)*(a(j+1,i+1)*z1-b(j+1,i+1)*z2);
            end
        end   
        if j==0
            s1(k)=delta11(j+1)*B/K*(-B/2+k*B/K-B/2/K);
        else
            z7=-(B/2/j/pi)^2*2*cos((2*j*pi)*(-1/2+k/K-1/2/K))*sin(2*j*pi/2/K);
            s1(k)=s1(k)+delta11(j+1)*z7;
        end
    end
    for i=0:N-1
        if i==0
            delta12(k)=delta12(k)-k2*c(1,i+1)*(phir(k)-phif(k));
        else
            z3=(sin(i*phir(k))-sin(i*phif(k)))/i;%中间变量
            z4=(cos(i*phir(k))-cos(i*phif(k)))/i;%中间变量
            delta12(k)=delta12(k)+(-i^2*k1(k)-k2)*(c(1,i+1)*z3-d(1,i+1)*z4);
        end
    end
    z8=(R+tau)*k1(k)*B/K*(cos(phif(k))-cos(phir(k)))+k2*B/K*jifen0(phif(k),phir(k),X(:,k));
    q1=q1+(0.5/pi/B)*(s1(k)+delta12(k)*B/K+z8);
end
%以上为求解路面激励的程序

%下面是求解带束层位移系数的程序        
z21=x(4);
z22=x(5);
z23=x(6);
z24=-c1/m1*x(4)-g/m1*x(5)-kmn/m1*x(1)-ka/m1*x(3)+q1/m1;
z25=-c1/m1*x(5)+g/m1*x(4)-kmn/m1*x(2);
z26=-ca/ma*x(6)-ka/ma*x(1)-kb/ma*x(3)+T/2/pi/R/B/ma;
z=[z21;z22;z23;z24;z25;z26];
end