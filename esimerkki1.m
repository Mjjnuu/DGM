close all;
clear all;

t=-pi:0.01:pi;
alpha=-1000;
beta=100;
u=2.*exp(-t).*(exp(pi)*sin(t./2) + cos(t./2)./exp(pi)*(alpha +1 -2*exp(2*pi))) + beta.*sin(t);
plot (t,u, 'LineWidth', 1.5);
hold on

%Yhtälö -d/dx(-4*du/dx) + 8*du/dx + 5*u = beta*(sinx + 8*cosx), -pi<x<pi
%Neumann: -4*du/dx = -4*alpha
%Dirichlet: u(pi) = 2

%Muodostetaan välille jako
tic
n=628;
x=zeros(1,n+1);
summa=0;
for i = 2:n+1
    summa = summa + rand;
    x(i)=x(i) + summa;
   
end
x = x./summa;
x=x-0.5;
x=x.*(2*pi);
A=zeros(n+1);
b=zeros(n+1,1);

for k = 1:n
    Alocal = zeros(2);
    blocal= zeros(2,1);
    h=x(k+1)-x(k);
    
    phi1 = @(q) (x(k+1)-q)/h;
    phi2 = @(q) (q-x(k))/h;
    
    %Kokeillaan ensimmäisen termin merkin vaihtamista
    a11 = @(q) -4/h.^2-8/h*phi1(q) + 5*phi1(q).^2; 
    a12 = @(q) 4/h.^2 +8/h*phi1(q) + 5*phi1(q).*phi2(q);
    a21 = @(q) 4/h.^2 -8/h*phi2(q) + 5*phi1(q).*phi2(q);
    a22 = @(q) -4/h.^2 +8/h*phi2(q) + 5*phi2(q).^2;
    
    b1 = @(q) beta*(sin(q) + 8*cos(q)).*phi1(q);
    b2 = @(q) beta*(sin(q) + 8*cos(q)).*phi2(q);
    
    Alocal(1,1)=integral(a11,x(k),x(k+1));
    Alocal(1,2)=integral(a12,x(k),x(k+1));
    Alocal(2,1)=integral(a21,x(k),x(k+1));
    Alocal(2,2)=integral(a22,x(k),x(k+1));
    blocal(1)=integral(b1,x(k),x(k+1));
    blocal(2)=integral(b2,x(k),x(k+1));
    
    A(k:k+1,k:k+1) = A(k:k+1,k:k+1) + Alocal;
    b(k:k+1) = b(k:k+1) + blocal;
    
    
end
A(n+1,:)=0;
A(n+1,n+1)=1;
b(n+1)=2;
b(1) = b(1) +4*alpha;

U=A\b;

disp(toc)
plot (x,U);
U(45)-u(45)



















