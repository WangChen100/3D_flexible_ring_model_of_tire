function [sy]=jifen1(phif,phir,X,n)
x=phif:(phir-phif)/50:phir;
y=(X(1)*x.^3+X(2)*x.^2+X(3)*x+X(4)).*cos(n*x);
sy=trapz(x,y);
end