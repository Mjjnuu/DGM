%Ratkaistaan yhtälö -d/dx(-x*du/dx)-(2x+3)du/dx+(x+2)u = x^3*e^(2x) ,
%0.5<x<1, reunaehtoina u(0.5)=e ja -xdu/dx=-5e^2 kun x=1.

%luodaan satunnaisesti vektori x joka on välin (0.5,1) jako
n=20;
x=zeros(1,n+1);
summa=0;
for i=2:n+1
    summa = summa + rand;
    x(i)=summa;
end
%normitetaan x, ja siirretään alkavaksi halutusta paikasta
x=x./(2*summa);
x=x+0.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t=0.5:0.01:1;
e=2.718281828459;
k1=3*e-32/31*(3*e+sqrt(e)/4);
k2=8/31*(3*e+sqrt(e)/4);
u=e.^t.*(k1+k2*t.^3)+e.^(2*t).*(t.^2-2*t+2);
plot (t,u)
hold on

tic
%Ratkaistaan finite element metodilla:
%Yhtälö muotoa -d/dx(p*du/dx)+q*du/dx+r*u=f
p = @(s) -s;
q = @(s) -(2*s+3);
r = @(s) s+2;
f = @(s) s.^3.*e.^(2*s);

%Alustetaan matriisi A ja vektori b
A=zeros(n+1);
b=zeros(n+1,1);

A_local=zeros(2);
b_local=zeros(2,1);

for k=1:n
    h=x(k+1)-x(k);
    phi1 = @(z) (x(k+1)-z)/(x(k+1)-x(k));
    phi2 = @(z) (z-x(k))/(x(k+1)-x(k));
    A11 = @(z) p(z)/h^2 - q(z).*phi1(z)/h + r(z).*phi1(z).^2;
    A12 = @(z) -p(z)/h^2 + q(z).*phi1(z)/h + r(z).*phi1(z).*phi2(z);
    A21 = @(z) -p(z)/h^2 - q(z).*phi2(z)/h + r(z).*phi2(z).*phi1(z);
    A22 = @(z) p(z)/h^2 + q(z).*phi2(z)/h + r(z).*phi2(z).^2;
    
    
    A_local(1,1)=integral(A11,x(k),x(k+1));
    A_local(1,2)=integral(A12,x(k),x(k+1));
    A_local(2,1)=integral(A21,x(k),x(k+1));
    A_local(2,2)=integral(A22,x(k),x(k+1));
    
    %Lisätään A:han alimatriisien arvot
    A(k:k+1,k:k+1)=A(k:k+1,k:k+1) + A_local;
    
    b1 = @(z) f(z).*phi1(z);
    b2 = @(z) f(z).*phi2(z);
    b_local(1)=integral(b1,x(k),x(k+1));
    b_local(2)=integral(b2,x(k),x(k+1));
    
    b(k:k+1)=b(k:k+1) + b_local;
end    
A(1,:)=0;
A(1,1)=1;
b(1)=e;
b(n+1)=b(n+1) + -5*exp(1)*exp(1); 

U=A\b;
disp(toc);
plot (x,U,'-o');
    
    
    
    
    
    
