a=-1;
b=5;
N=2000;
t=a:(b-a)/N:b;

%The number of points on a partition of the domain
n=16;
%Partition x
x=zeros(1,n);
Sum = 0;
for i =2:n
    Sum = Sum + rand + 0.25;
    x(i)=Sum;
end
x=x./Sum.*(b-a)+a;


y=zeros(1,N+1);
tic
for i=1:n-1
    [A,l,u]=GetDomain(t,x(i),x(i+1));
    y(l:u)=LegendreTransform(i,A)./sqrt(i); 
   
end    
plot(t,y);
disp(toc);
hold on