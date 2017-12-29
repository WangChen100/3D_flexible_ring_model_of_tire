function z=myfun2_2(t,x)
%计算m!=0,n!=0是的系数a、b
global ringdensity omega R B h p0 sigmah sigmahy kw kv ku lamida miu Eh D G k1 k2 k3 m n M N K a b phif phir;%定义全局变量；
m1=ringdensity*(n^2+2);
g=2*n^3*ringdensity*omega;
kmn=ringdensity*omega^2*(-n^4+n^2-1)+sigmah/R^2*(n^4+n^2+1)-p0/R+(D/R^4)*(n^6-2*n^4+n^2)+n^2*kw+kv+ku-...
    (2*miu*D/R^2)*(n^4-n^2)*(2*m*pi/B)^2+D*n^2*(2*m*pi/B)^4+sigmahy*(n^2+1)*(2*m*pi/B)^2+G*h^3/12/R^3*(4*n^4-4*n^2+1)*(2*m*pi/B)^2+Eh/(1-miu^2)*(2*m*pi/B)^2;
c1=2*lamida*sqrt(m1*kmn);


delta2=zeros(1,M);
delta3=zeros(1,M);
s2=zeros(1,K);
s3=zeros(1,K);
q2=0;%n=0路面合激励
q3=0;%n!=0xsinn路面合激励
for k=1:K
    for j=0:M-1
        for i=0:N-1
            if i==n
                    z1=(phir(k)-phif(k)+eps)/2+(sin((i+n)*phir(k))-sin((i+n)*phif(k))+eps)/2/(i+n);%中间变量
                    z2=-(cos((i+n)*phir(k))-cos((i+n)*phif(k))+eps)/2/(i+n);%中间变量
                    z3=-(cos((i+n)*phir(k))-cos((i+n)*phif(k))+eps)/2/(i+n);%中间变量
                    z4=(phir(k)-phif(k)+eps)/2-(sin((i+n)*phir(k))-sin((i+n)*phif(k))+eps)/2/(i+n);%中间变量
                    delta2(j+1)=delta2(j+1)+(-i^2*k1(k)-k2-k3)*(a(j+1,i+1)*z1+b(j+1,i+1)*z2);   
                    delta3(j+1)=delta3(j+1)+(-i^2*k1(k)-k2-k3)*(a(j+1,i+1)*z3+b(j+1,i+1)*z4);
            else
                    z5=(sin((i-n)*phir(k))-sin((i-n)*phif(k))+eps)/2/(i-n)+(sin((i+n)*phir(k))-sin((i+n)*phif(k))+eps)/2/(i+n);%中间变量
                    z6=-(cos((i+n)*phir(k))-cos((i+n)*phif(k))+eps)/2/(i+n)-(cos((i-n)*phir(k))-cos((i-n)*phif(k))+eps)/2/(i-n);%中间变量        
                    z7=-(cos((i+n)*phir(k))-cos((i+n)*phif(k))+eps)/2/(i+n)+(cos((i-n)*phir(k))-cos((i-n)*phif(k))+eps)/2/(i-n);%中间变量
                    z8=(sin((i-n)*phir(k))-sin((i-n)*phif(k))+eps)/2/(i-n)-(sin((i+n)*phir(k))-sin((i+n)*phif(k))+eps)/2/(i+n);%中间变量
                    delta2(j+1)=delta2(j+1)+(-i^2*k1(k)-k2-k3)*(a(j+1,i+1)*z5+b(j+1,i+1)*z6); 
                    delta3(j+1)=delta3(j+1)+(-i^2*k1(k)-k2-k3)*(a(j+1,i+1)*z7+b(j+1,i+1)*z8);
            end
        end
        if j==0
            s2(k)=delta2(j+1)*(-B/2/m/pi)*(cos(2*m*pi*(-1/2+k/K))-cos(2*m*pi*(-1/2+k/K-1/K)));
            s3(k)=delta3(j+1)*(-B/2/m/pi)*(cos(2*m*pi*(-1/2+k/K))-cos(2*m*pi*(-1/2+k/K-1/K)));
        else
            if j==m
                z9=B/2/K-B/2/(j+m)/pi*sin((j+m)*pi/K)*cos(2*(j+m)*pi*(-1/2+k/K-1/2/K));
                s2(k)=s2(k)+delta2(j+1)*z9;
                s3(k)=s3(k)+delta3(j+1)*z9;
            else
                z10=B/2/(j-m)/pi*sin((j-m)*pi/K)*cos(2*(j-m)*pi*(-1/2+k/K-1/2/K))-B/2/(j+m)/pi*sin((j+m)*pi/K)*cos(2*(j+m)*pi*(-1/2+k/K-1/2/K));
                s2(k)=s2(k)+delta2(j+1)*z10;
                s3(k)=s3(k)+delta3(j+1)*z10;
            end
        end
    end
    q2=q2+(2/pi/B)*s2(k);
    q3=q3+(2/pi/B)*s3(k);
end
%以上为求解路面激励的程序

%下面是求解带束层位移系数的程序        
z21=x(3);
z22=x(4);
z23=-c1/m1*x(3)-g/m1*x(4)-kmn/m1*x(1)+q2/m1;
z24=-c1/m1*x(4)+g/m1*x(3)-kmn/m1*x(2)+q3/m1;
z=[z21;z22;z23;z24];
end