function [z1 z2 z3 z4]=bjtj(phi,fx)
global L;
z1=0;
z2=0;
if fx((L+1)/2)<=0
    z1=0;
    z2=0;
    z3=(L+1)/2;
    z4=(L+1)/2;
    return;
end
if fx(1)>0||fx(end)>0
    disp(['error!stress is too high']);
    return;
end
for i=1:(L-1)/2
    if fx((L+1)/2-i)<=0
        z1=phi((L+1)/2-i);
        z3=(L+1)/2-i;
        break;
    end
end
for j=1:(L-1)/2
    if fx((L+1)/2+j)<=0
        z2=phi((L+1)/2+j);
        z4=(L+1)/2+j;
        break;
    end
end
end
 