function [X]=vnihe(phif,phir,phim,v1,v2,v3)
if phif~=phir
    if v1==v3
        AA=[phif^3,phif^2,phif,1;phir^3,phir^2,phir,1;
            3*phir^2,2*phir,1,0;3*phif^2,2*phif,1,0];%需要修改
        BB=[v1;v2;-1e3;0];
        X=AA\BB;
    elseif v2==v3
        AA=[phif^3,phif^2,phif,1;phir^3,phir^2,phir,1;
            3*phif^2,2*phif,1,0;3*phir^2,2*phir,1,0];%需要修改
        BB=[v1;v2;-1e3;0];
        X=AA\BB;
    else
        AA=[phif^3,phif^2,phif,1;phir^3,phir^2,phir,1;
            phim^3,phim^2,phim,1;3*phim^2,2*phim,1,0];
        BB=[v1;v2;v3;0];
        X=AA\BB;
    end
else
    X=0;
end
