clc
clear all
d=0;vd=0;ztx=0;vztx=0;
[t f]=ode23(@odefun3,0:0.01:0.4,[d;ztx;vd;vztx]);
% ff(1:11,:)=f;
d=f(end,1);
ztx=f(end,2);
vd=f(end,3);
vztx=f(end,4);