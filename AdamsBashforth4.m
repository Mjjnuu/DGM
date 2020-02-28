function results=AdamsBashforth4(f,t,RESULTS)
    n=length(t)-1;
    
    Memory=zeros(size(RESULTS,1),4);
    results=RESULTS;
    results(:,1:4)=RungeKutta4(f,t(1:4),RESULTS(:,1:4));
    Memory(:,2)=(t(2)-t(1))*f(t(1),results(1));
    Memory(:,3)=(t(3)-t(2))*f(t(2),results(2));
    Memory(:,4)=(t(4)-t(3))*f(t(3),results(3));
    
    for i=4:n
        Memory(:,1)=Memory(:,2);
        Memory(:,2)=Memory(:,3);
        Memory(:,3)=Memory(:,4);
        Memory(:,4)=(t(i+1)-t(i))*f(t(i),results(:,i));
        results(i+1)=results(i)+55/24*Memory(:,4)-59/24*Memory(:,3)+37/24*Memory(:,2)-9/24*Memory(:,1);
        
    end    
end