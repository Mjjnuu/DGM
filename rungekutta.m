close all;
clear all;

a=0;
b=10;
n=100;
h=(b-a)/n;
t=a:h:b;
%yhtälö y'=y-3*sin(3*t)/(t+1) alkuehto y(0)=1

%Eulerin metodi
y=ones(1,n+1);
tic
for i = 1:n
    %y(t(i+1))=y(t(i))+y'(t(i))*h
    y(i+1)=y(i)+h*f(t(i),y(i));
    
end
disp(" ");
disp("Eulerin metodi: "+ toc + ", red");
plot (t,y,'r');%%%%%%%%%%%%%%%%%%%%%%%%%

hold on

%Midpoint-metodi
%y(n+1)=y(n)+hf( t(n)+h/2, y(n) + h/2*f(t(n),y(n))
Y=ones(1,n+1);
tic
for i= 1:n
    Y(i+1)=Y(i)+h*f( t(i)+h/2 , Y(i)+h/2*f(t(i),Y(i)));
end
disp("Midpoint: " + toc + ", blue")
%plot (t,Y,'b');%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%


%Runge-kutta RK4
%y(n+1)=y(n)+1/6*(k1+2*k2+2*k3+k4)
%k1 = h*f(t(n), y(n))
%k2 = h*f(t(n)+h/2, y(n)+k1/2)
%k3 = h*f(t(n)+h/2, y(n)+k2/2)
%k4 = h*f(t(n)+h, y(n)+k3)
z=ones(1,n+1);
tic
for i = 1:n
    k1=h*f(t(i), z(i));
    pause(0.001)
    k2=h*f(t(i)+h/2, z(i)+k1/2);
    pause(0.001)
    k3=h*f(t(i)+h/2, z(i)+k2/2);
    pause(0.001)
    k4=h*f(t(i)+h,z(i)+k3);
    pause(0.001)
    z(i+1)=z(i)+1/6*(k1+2*k2+2*k3+k4);
end
disp("Runge-Kutta: " + toc + ", dashed")
plot (t,z,'--');

%Adams-Bashfort method
w=ones(1,n+1);
muisti=zeros(1,2);
tic

k1=h*f(t(1), w(1));
k2=h*f(t(1)+h/2, w(1)+k1/2);
k3=h*f(t(1)+h/2, w(1)+k2/2);
k4=h*f(t(1)+h,w(1)+k3);
w(2)=w(1)+1/6*(k1+2*k2+2*k3+k4);

muisti(2)=h*f(t(1),w(1));

for i = 2:n
    muisti(1)=muisti(2);
    muisti(2)=h*f(t(i),w(i));
    %w(i+1) = w(i) + 3/2*h*f(t(i),w(i))-1/2*h*f(t(i-1),w(i-1));
    w(i+1)=w(i) + 3/2*muisti(2)-1/2*muisti(1);
end    

disp("Adams-Bashforth: " + toc + ", green")
%plot (t,w,'g');%%%%%%%%%%%%%%%%%%%5



results=ones(1,n/2);
results=RungeKutta4(@f,t(1:2:end),results);
plot(t(1:2:end),results,'-o' )



%4th order Adams-Basforth method
%y(n+4)=y(n+3)+h*(55/24*f(n+3)-59/24*f(n+2)+37/24*f(n+1)-9/24*f(n))
u=ones(1,n+1);
Muisti=zeros(1,4);
tic
for i=1:3
    k1=h*f(t(i), u(i));
    pause(0.001)
    k2=h*f(t(i)+h/2, u(i)+k1/2);
    pause(0.001)
    k3=h*f(t(i)+h/2, u(i)+k2/2);
    pause(0.001)
    k4=h*f(t(i)+h,u(i)+k3);
    pause(0.001)
    u(i+1)=u(i)+1/6*(k1+2*k2+2*k3+k4);
end
Muisti(2)=h*f(t(1),u(1));
pause(0.001)
Muisti(3)=h*f(t(2),u(2));
pause(0.001)
Muisti(4)=h*f(t(3),u(3));
pause(0.001)
for i=4:n
    Muisti(1)=Muisti(2);
    Muisti(2)=Muisti(3);
    Muisti(3)=Muisti(4);
    Muisti(4)=h*f(t(i),u(i));
    pause(0.001)
    u(i+1)=u(i)+55/24*Muisti(4)-59/24*Muisti(3)+37/24*Muisti(2)-9/24*Muisti(1);
end
disp("4th order A-B: " + toc + ", dashed red")
plot (t,u,'--r');%%%%%%%%%%%%%%%%%%%







disp(" ")




N=100000;
h=(b-a)/N;
t=a:h:b;
%t(50001)
%yhtälö y'=y-3*sin(3*t)/(t+1) alkuehto y(0)=1

%Eulerin metodi
y=ones(1,N+1);
tic
for i = 1:N
    k1=h*f(t(i), y(i));
    k2=h*f(t(i)+h/2, y(i)+k1/2);
    k3=h*f(t(i)+h/2, y(i)+k2/2);
    k4=h*f(t(i)+h,y(i)+k3);
    y(i+1)=y(i)+1/6*(k1+2*k2+2*k3+k4);
end
%fprintf('%.32f', y(50001));
%disp(toc)
plot (t,y,'k');




    
function z=f(t,y)
    %z=y^2-y-3*sin(3*t)/(t+1);%alkuperäinen
    z=y^2-y-3*sin(t*t);
    %z=sin((2*pi)*t)
    
end
