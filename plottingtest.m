clear all;
a=-2;
b=2;
d=b-a;
resolution=100;
t=a:1/resolution:b;
tic
N=10;

S=0;
x=zeros(1,N+1);
for i = 2:N+1
    S=S+rand+0.2;
    x(i)=S;
end
x=x./S.*d+a;

%fun = @(q) sin(4*q)+log(q.^2+2)-0.001*q.^10 + sin(q+2*cos(q+3*sin(q+4*cos(q))));
fun = @(q) exp(q)/3-sin(2.*q)

plot(t,fun(t))
hold on
Z=ApprLP(fun,t,x,4);
toc
%Y=plotLP(t,x,ApprLP(fun,t,x,2));
tic
Y=plotLP(t,x,Z);

plot(t,Y)
disp(" ")
disp(toc + "s")
disp(" ")
