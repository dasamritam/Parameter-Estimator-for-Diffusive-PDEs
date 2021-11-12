function [value] = BBspline(d,i,t,T)
    if(d==0)
        if(t>T(i) && t<=T(i+1))
            value=1;
            return;
        else
            value=0;
            return;
        end
    end
    value=(t-T(i))/(T(i+d)-T(i))*BBspline(d-1,i,t,T)+...
        (T(i+d+1)-t)/(T(i+d+1)-T(i+1))*BBspline(d-1,i+1,t,T);
    if(isnan(value))
            value=0;
    end
    return;
end
 