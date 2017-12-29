function z=odefun3(t,x)
m1=22.4;
k1=1.1326e6;
ma=-1.5417;
k2=9.2605e4;
ka=-2.9873e5;
lamida=0;%×èÄá±È
c1=2*lamida*sqrt(k1/m1);
ca=2*lamida*sqrt(k2/ma);
z1=x(3);
z2=x(4);
z3=-c1*x(3)-k1/m1*x(1)-ka/m1*x(2)+0/m1;%
z4=-ca*x(4)-ka/ma*x(1)-k2/ma*x(2)+0.0023;
z=[z1;z2;z3;z4];
end