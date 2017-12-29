function [sy]=jifen0(phif,phir,X)
a=phir^4-phif^4;
b=phir^3-phif^3;
c=phir^2-phif^2;
d=phir-phif;
sy=X(1)/4*a+X(2)/3*b+X(3)/2*c+X(4)*d;
end