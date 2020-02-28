clear all;
a=-2;
b=4;
L=b-a;
resolution=50;
t=a:1/resolution:b;

%fn = @(w) exp(w./2)-cos(2*w)+sin(5*w)./2-w;
%fn = @(w) sin(1./(w.^2+0.3));
%fn = @(w) (exp(1./(.5*w.^4+1)) /2).^2;
%fn = @(w) sqrt(0.5-(-.2-w).^2)./3;
%fn = @(w) ( ( exp( 1./((-0.2-w).^2+1) ) ).^4  )./100;
fn = @(w) (exp(1./(.5*w.^4+1)) /2).^2 + ...
    ( ( exp( 1./((-0.15-w).^2+1) ) ).^4  )./150 + ...
    exp(-(2.*w-3.8).^2)./3 + 0.02*(w-0.9).^2;
%fn = @(w) exp(-(2.*w-4).^2)./3;
%fn = @(w) 0.02*(w+1).^2;


plot(t,fn(t))
hold on

%The number of elements
N=5;
%The highest degree of Legendre polynomials
P=6;


Sum=0;
x=zeros(1,N+1);
for i = 2:N+1
    Sum=Sum+rand+0.2;
    x(i)=Sum;
end
x=x./Sum.*L+a;


LP=cell(26,1);
for i=1:26
    LP{i} = @(w) legendreP(i-1,w);
end


LPfast=cell(10,1);
LPfast{1} = @(x) x;
LPfast{2} = @(x) 1/2*(3*x.^2-1);
LPfast{3} = @(x) 1/2*(5*x.^3-3.*x);
LPfast{4} = @(x) 1/8*(35*x.^4-30*x.^2+3);
LPfast{5} = @(x) 1/8*(63*x.^5-70*x.^3+15*x);
LPfast{6} = @(x) 1/16*(231*x.^6-315*x.^4+105*x.^2-5);
LPfast{7} = @(x) 1/16*(429*x.^7-693*x.^5+315*x.^3-35*x);
LPfast{8} = @(x) 1/128*(6435*x.^8-12012*x.^6+6930*x.^4-1260*x.^2+35);
LPfast{9} = @(x) 1/128*(12155*x.^9-25740*x.^7+18018*x.^5-4620*x.^3+315*x);
LPfast{10} = @(x) 1/256*(46189*x.^10-109395*x.^8+90090*x.^6-30030*x.^4+3465*x.^2-63);


diagonal=2./(2*(0:P)+1);
Diagonal=zeros(1,(P+1)*N);

%(P+1)xN matrix containing integrals F_{ij}=\int_j f*L_i dx
F=zeros(P+1,N);



for i=1:N
    %deltaxi=x(i+1)-x(i);
    for power = 0:P
        
        F(power + 1, i)=LegendreProduct(fn,power,x(i),x(i+1));
        
    end 
    Diagonal((i-1)*(P+1)+1  : i*(P+1))=diagonal;%deltaxi*diagonal;
end

F=reshape(F,[],1);


%Mass matrix
I=(1:(P+1)*N);
M=sparse(I, I, Diagonal);
%M=diag(Diagonal);
U=M\F;
U=reshape(U,P+1,N);




%Plotting%%%%%%%%%%%%%%

Y=zeros(1,length(t));
for element = 1:N
    %
    tic
    [subdomain,l,u] = GetDomain2(t,x(element),x(element+1));
    values=zeros(1,length(subdomain));
    for p=0:P
       values=values + LegendreTransform(p,subdomain)*U(p+1,element);
       
    end    
    Y(l:u)=values;
    disp(toc)
end    
%disp("Plotting time: " + toc)
plot(t,Y)

disp("L2 error: " + sum(((fn(t)-Y)./resolution).^2))















