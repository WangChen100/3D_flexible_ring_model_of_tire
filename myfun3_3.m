function z=myfun3_3(t,x)
%计算系数c、d
global ringdensity omega R B p0 sigmah kw kv lamida tau D k1 k2 n M N K a b c d phif phir v1 v2 v3 X;%定义全局变量；
m1=ringdensity*(n^2+1);
g=2*(n^3-n)*ringdensity*omega;
kmn=-ringdensity*omega^2*(n^2-1)^2+p0/R*(n^2-1)+(D/R^4)*(n^6-2*n^4+n^2)+(sigmah/R^2)*(n^4+1)+n^2*kw+kv;
c1=2*lamida*sqrt(m1*kmn);


delta21=zeros(1,M);
delta22=zeros(1,K);
delta23=zeros(1,K);
delta31=zeros(1,M);
delta32=zeros(1,K);
delta33=zeros(1,K);

s2=zeros(1,K);
s3=zeros(1,K);

q2=0;%n!=0xcosn路面合激励
q3=0;%n!=0xsinn路面合激励
for k=1:K
    delta23(k)=-(cos((n+1)*phir(k))-cos((n+1)*phif(k)))/(n+1)+(cos((n-1)*phir(k))-cos((n-1)*phif(k)))/(n-1);
    delta33(k)=(sin((1-n)*phir(k))-sin((1-n)*phif(k)))/(1-n)-(sin((n+1)*phir(k))-sin((n+1)*phif(k)))/(n+1);

    for j=0:M-1
        for i=0:N-1
                if i==n
                    z1=(phir(k)-phif(k))/2+(sin((i+n)*phir(k))-sin((i+n)*phif(k)))/2/(i+n);%中间变量
                    z2=-(cos((i+n)*phir(k))-cos((i+n)*phif(k)))/2/(i+n);%中间变量
                    z5=-(cos((i+n)*phir(k))-cos((i+n)*phif(k)))/2/(i+n);%中间变量
                    z6=(phir(k)-phif(k))/2-(sin((i+n)*phir(k))-sin((i+n)*phif(k)))/2/(i+n);%中间变量
                    delta21(j+1)=delta21(j+1)+(-i^2*k1(k)-k2)*(a(j+1,i+1)*z1+b(j+1,i+1)*z2);
                    delta31(j+1)=delta31(j+1)+(-i^2*k1(k)-k2)*(a(j+1,i+1)*z5+b(j+1,i+1)*z6);
                else
                    z3=(sin((i-n)*phir(k))-sin((i-n)*phif(k)))/2/(i-n)+(sin((i+n)*phir(k))-sin((i+n)*phif(k)))/2/(i+n);%中间变量
                    z4=-(cos((i+n)*phir(k))-cos((i+n)*phif(k)))/2/(i+n)-(cos((i-n)*phir(k))-cos((i-n)*phif(k)))/2/(i-n);%中间变量
                    z7=-(cos((i+n)*phir(k))-cos((i+n)*phif(k)))/2/(i+n)+(cos((i-n)*phir(k))-cos((i-n)*phif(k)))/2/(i-n);%中间变量
                    z8=(sin((i-n)*phir(k))-sin((i-n)*phif(k)))/2/(i-n)-(sin((i+n)*phir(k))-sin((i+n)*phif(k)))/2/(i+n);%中间变量
                    delta21(j+1)=delta21(j+1)+(-i^2*k1(k)-k2)*(a(j+1,i+1)*z3+b(j+1,i+1)*z4);
                    delta31(j+1)=delta31(j+1)+(-i^2*k1(k)-k2)*(a(j+1,i+1)*z7+b(j+1,i+1)*z8);
                end
        end   
        if j==0
            s2(k)=delta21(j+1)*B/K*(-B/2+k*B/K-B/2/K);
            s3(k)=delta31(j+1)*B/K*(-B/2+k*B/K-B/2/K);
        else
            z9=-(B/2/j/pi)^2*2*cos((2*j*pi)*(-1/2+k/K-1/2/K))*sin((2*j*pi)*(1/2/K));
            s2(k)=s2(k)+delta21(j+1)*z9;
            s3(k)=s3(k)+delta31(j+1)*z9;
        end
    end
    for i=0:N-1
        if i==n
            z1=(phir(k)-phif(k))/2+(sin((i+n)*phir(k))-sin((i+n)*phif(k)))/2/(i+n);%中间变量
            z2=-(cos((i+n)*phir(k))-cos((i+n)*phif(k)))/2/(i+n);%中间变量
            z5=-(cos((i+n)*phir(k))-cos((i+n)*phif(k)))/2/(i+n);%中间变量
            z6=(phir(k)-phif(k))/2-(sin((i+n)*phir(k))-sin((i+n)*phif(k)))/2/(i+n);%中间变量
            delta22(k)=delta22(k)+(-i^2*k1(k)-k2)*(c(1,i+1)*z1+d(1,i+1)*z2);
            delta32(k)=delta32(k)+(-i^2*k1(k)-k2)*(c(1,i+1)*z5+d(1,i+1)*z6);
        else
            z3=(sin((i-n)*phir(k))-sin((i-n)*phif(k)))/2/(i-n)+(sin((i+n)*phir(k))-sin((i+n)*phif(k)))/2/(i+n);%中间变量
            z4=-(cos((i+n)*phir(k))-cos((i+n)*phif(k)))/2/(i+n)-(cos((i-n)*phir(k))-cos((i-n)*phif(k)))/2/(i-n);%中间变量
            z7=-(cos((i+n)*phir(k))-cos((i+n)*phif(k)))/2/(i+n)+(cos((i-n)*phir(k))-cos((i-n)*phif(k)))/2/(i-n);%中间变量
            z8=(sin((i-n)*phir(k))-sin((i-n)*phif(k)))/2/(i-n)-(sin((i+n)*phir(k))-sin((i+n)*phif(k)))/2/(i+n);%中间变量
            delta22(k)=delta22(k)+(-i^2*k1(k)-k2)*(c(1,i+1)*z3+d(1,i+1)*z4);
            delta32(k)=delta32(k)+(-i^2*k1(k)-k2)*(c(1,i+1)*z7+d(1,i+1)*z8);
        end
    end
%     z9=(v2(k)*sin(n*phir(k))-v1(k)*sin(n*phif(k)))/n;
%     z10=v3(k)*(cos(n*phir(k))-cos(n*phif(k)))/n^2;
%     z11=(v1(k)*cos(n*phif(k))-v2(k)*cos(n*phir(k)))/n;
%     z12=v3(k)*(sin(n*phir(k))-sin(n*phif(k)))/n^2;
    z9=jifen1(phif(k),phir(k),X(:,k),n);
    z11=jifen2(phif(k),phir(k),X(:,k),n);
    z13=(R+tau)*k1(k)*B/K/2*delta23(k)+k2*(z9)*B/K;%z9+z10
    z14=(R+tau)*k1(k)*B/K/2*delta33(k)+k2*(z11)*B/K;%z11+z12
    q2=q2+(1/pi/B)*(s2(k)+delta22(k)*B/K+z13);
    q3=q3+(1/pi/B)*(s3(k)+delta32(k)*B/K+z14);
end
%以上为求解路面激励的程序

%下面是求解带束层位移系数的程序        
z21=x(3);
z22=x(4);
z23=-c1/m1*x(3)-g/m1*x(4)-kmn/m1*x(1)+q2/m1;
z24=-c1/m1*x(4)+g/m1*x(3)-kmn/m1*x(2)+q3/m1;
z=[z21;z22;z23;z24];
end