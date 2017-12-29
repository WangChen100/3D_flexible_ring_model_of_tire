function z=myfun2_1(t,x)
%����m!=0,n=0�ǵ�ϵ��a��b
global ringdensity omega R B h p0 sigmah sigmahy kw kv ku lamida miu Eh D G k1 k2 k3 m n M N K a b phif phir;%����ȫ�ֱ�����
m1=ringdensity*(n^2+2);
g=2*n^3*ringdensity*omega;
kmn=ringdensity*omega^2*(-n^4+n^2-1)+sigmah/R^2*(n^4+n^2+1)-p0/R+(D/R^4)*(n^6-2*n^4+n^2)+n^2*kw+kv+ku-...
    (2*miu*D/R^2)*(n^4-n^2)*(2*m*pi/B)^2+D*n^2*(2*m*pi/B)^4+sigmahy*(n^2+1)*(2*m*pi/B)^2+G*h^3/12/R^3*(4*n^4-4*n^2+1)*(2*m*pi/B)^2+Eh/(1-miu^2)*(2*m*pi/B)^2;
c1=2*lamida*sqrt(m1*kmn);


delta1=zeros(1,M);
s1=zeros(1,K);
q1=0;%n=0·��ϼ���
for k=1:K
    for j=0:M-1
        for i=0:N-1
            if i==0
                  delta1(j+1)=delta1(j+1)+(-k2-k3)*a(j+1,i+1)*(phir(k)-phif(k));
            else
                  z1=sin(i*phir(k))-sin(i*phif(k));%�м����
                  z2=cos(i*phir(k))-cos(i*phif(k));%�м����
                  delta1(j+1)=delta1(j+1)+(-i^2*k1(k)-k2-k3)/i*(a(j+1,i+1)*z1-b(j+1,i+1)*z2);
            end
        end
        if j==0
            s1(k)=delta1(j+1)*(-B/2/m/pi)*(cos(2*m*pi*(-1/2+k/K))-cos(2*m*pi*(-1/2+k/K-1/K)));
        else
            if j==m
                z7=B/2/K-B/2/(j+m)/pi*sin((j+m)*pi/K)*cos(2*(j+m)*pi*(-1/2+k/K-1/2/K));
                s1(k)=s1(k)+delta1(j+1)*z7;
            else
                z8=B/2/(j-m)/pi*sin((j-m)*pi/K)*cos(2*(j-m)*pi*(-1/2+k/K-1/2/K))-B/2/(j+m)/pi*sin((j+m)*pi/K)*cos(2*(j+m)*pi*(-1/2+k/K-1/2/K));
                s1(k)=s1(k)+delta1(j+1)*z8;
            end
        end
    end
    q1=q1+(1/pi/B)*s1(k);
end
%����Ϊ���·�漤���ĳ���

%��������������λ��ϵ���ĳ���        
z21=x(3);
z22=x(4);
z23=-c1/m1*x(3)-g/m1*x(4)-kmn/m1*x(1)+q1/m1;
z24=-c1/m1*x(4)+g/m1*x(3)-kmn/m1*x(2);
z=[z21;z22;z23;z24];
end