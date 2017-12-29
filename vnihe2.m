function [z4]=vnihe2(phi,phif,phir,phim,v1,v2,v3,z1,z2,z3)
if phif~=phir
    if z1==z3
        kb=(v2-v3)/(phir-phim)^2;
        z4(z3-z1+2:z2-z1+1)=kb*(phi(z3+1:z2)-phim).^2+v3;
    elseif z2==z3
        ka=(v1-v3)/(phif-phim)^2;
        z4(1:z3-z1+1)=ka*(phi(z1:z3)-phim).^2+v3;
    else
        ka=(v1-v3)/(phif-phim)^2;
        z4(1:z3-z1+1)=ka*(phi(z1:z3)-phim).^2+v3;
        kb=(v2-v3)/(phir-phim)^2;
        z4(z3-z1+2:z2-z1+1)=kb*(phi(z3+1:z2)-phim).^2+v3;
    end
else
    z4=0;
end
