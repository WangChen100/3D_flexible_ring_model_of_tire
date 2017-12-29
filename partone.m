function [XX]=partone(R,tau,N,Phir,Phif,rn,An,kk2)
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
