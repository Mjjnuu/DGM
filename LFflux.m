function y=LFflux(F,dF,a,b)
    %C=max(s \in [min(a,b),max(a,b)], |F'(s)|)
    %C is the maximum value of |F'(s)| between a and b.
    %vector containing all the lower bounds
    l=min(a,b);%unnecessary
    %vector containing all the upper bounds
    u=max(a,b);%unnecessary
    %C=max(abs(dF(l:(u-l)/10:u)));
    C=max(abs(dF(l+kron((0:10),(u-l)./10))), [], 2);
    y=1/2*(F(a)+F(b)-C.*(b-a));


end