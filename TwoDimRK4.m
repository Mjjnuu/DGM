function [resultsA,resultsU]=TwoDimRK4(f,g,t,RESULTSA,RESULTSU)
    n=length(t)-1;
    resultsA=RESULTSA;
    resultsU=RESULTSU;
    for i=1:n
        T=t(i);
        x=resultsA(:,i);
        y=resultsU(:,i);
        h=t(i+1)-t(i);
        
        k1a=h*f(T,x,y);
        k1u=h*g(T,x,y);
        
        k2a=h*f(T+h/2,x+k1a/2,y+k1u/2);
        k2u=h*g(T+h/2,x+k1a/2,y+k1u/2);
        
        k3a=h*f(T+h/2,x+k2a/2,y+k2u/2);
        k3u=h*g(T+h/2,x+k2a/2,y+k2u/2);
        
        k4a=h*f(T+h,x+k3a,y+k3u);
        k4u=h*g(T+h,x+k3a,y+k3u);
        
        resultsA(:,i+1)= x + 1/6*(k1a+2*k2a+2*k3a+k4a);
        resultsU(:,i+1)= y + 1/6*(k1u+2*k2u+2*k3u+k4u);
    end    
end