%Yhtälönä -d^2u/dx^2=2, 0<x<1, u(0)=u(1)=0
close all;
clear all;

t=0:0.01:1;
n=10;
h=1/n;
x=0:h:1;
A=zeros(n+1);
b=zeros(n+1,1);

for k=1:n
    Alocal=[1/h,-1/h; -1/h,1/h];
    blocal=[h;h];
    
    A(k:k+1,k:k+1) = A(k:k+1, k:k+1) + Alocal;
    b(k:k+1) = b(k:k+1) + blocal;
end

%Asetetaan reuna-arvot:
A(1,:)=0;%Ensimmäinen rivi nollia, lukuunottamatta indeksiä (1,1)
A(1,1)=1;
A(n+1,:)=0;
A(n+1,n+1)=1;
b(1)=0;
b(n+1)=0;

disp(A);
%Ratkaistaan kertoimet U:
%AU=b ==> U=A\b
U=A\b;

%Koska kantafunktiot ovat lineaarisia, niiden summa on lineaarinen.
%Matlabin plot piirtää pisteiden väliin lineaarisen interpolaation,
%joten ei tarvitse käsitellä kantafunktioita eksplisiittisesti:
plot (x,U);

hold on
y=t-t.^2;
plot (t,y);


