function [modms]=modifyms(meanms,y)%×İ²àÏòÈıÎ¬Ä¥Ëğ
mf=zeros(1,15);
modms=zeros(1,15);
for i=1:15
    if y(i)<-0.04
        mf(i)=1+(1-0.8)/(0.0825-0.04)*(y(i)+0.04);
    elseif y(i)>0.04
        mf(i)=1-(1-0.8)/(0.0825-0.04)*(y(i)-0.04);
    else
        mf(i)=1;
    end
    modms(i)=mf(i)*meanms(i);
end
end