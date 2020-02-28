function [t,lower, upper] = GetDomain(T,a,b)
    %T input vector, components assumed to be uniformally spaced
    %a lower boundary
    %b upper boundary
    %Returns t \subseteq T, which satisfies for all i=1,...,length(T)
    %a<=t(i)<=b
    
    n=length(T);
    lower=1;
    upper=n;
    if a<b
        if T(1)<=a
            if a<T(n)
                %a is checked
                %
                if b<=T(n)
                    if T(1)<b
                        %
                        %
                        %Both boundaries are lower than the real ones,
                        %but if lower bound is exact, it is used, unlike
                        %the upper bound. 
                        lower=floor( (n-1)*(a-T(1))/(T(n)-T(1)) ) + 1 ;
                        upper=floor( (n-1)*(b-T(1))/(T(n)-T(1)) ) + 1 ;
                        %
                        
                        t=T(lower:upper);
                        %toc
                        %
                        %
                    else
                        disp("ERROR: b is too small");
                    end
        
                else 
                    disp("ERROR: b is too large");
                end    
                %
                %
            else
                disp("ERROR: a is too large");
            end
        
        else 
            disp("ERROR: a is too small");
        end    
    else
        disp("Error: a is supposed to be smaller than b.");
    end
    
   

end