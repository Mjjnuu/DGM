function [t,lower, upper] = GetDomain2(T,a,b)
    
    n=length(T);
    lower=floor( (n-1)*(a-T(1))/(T(n)-T(1)) + 1 ) ;
    upper=floor( (n-1)*(b-T(1))/(T(n)-T(1)) + 1 ) ;
    %toc
    t=T(lower:upper);
    
end