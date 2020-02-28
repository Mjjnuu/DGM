function results=RungeKutta4(f,t,RESULTS)
    n=length(t)-1;
    results=RESULTS;
    for i=1:n
        T=t(i);
        x=results(:,i);
        h=t(i+1)-t(i);
        k1=h*f(T,x);
        k2=h*f(T+h/2,x+k1/2);
        k3=h*f(T+h/2,x+k2/2);
        k4=h*f(T+h,x+k3);
        
        results(:,i+1)= x + 1/6*(k1+2*k2+2*k3+k4);
    end    
end