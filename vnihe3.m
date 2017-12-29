function [z4]=vnihe3(phi,phif,phir,phim,v1,v2,v3)
if phif~=phir
    if phif==phim
        x=[phim (phim+phir)/2 phir];
        y=[v3 (v3+v2)/2*0.2 v2];
        z4=interp1(x,y,phi,'cubic');%spline
    elseif phir==phim
        x=[phif (phif+phim)/2 phim];
        y=[v1 (v1+v3)/2*0.2 v3];
        z4=interp1(x,y,phi,'cubic');
    else
        x=[phif (phif+phim)/2 phim (phim+phir)/2 phir];
        y=[v1 (v1+v3)/2*0.2 v3 (v3+v2)/2*0.2 v2];
        z4=interp1(x,y,phi,'cubic');
    end
else
    z4=0;
end
