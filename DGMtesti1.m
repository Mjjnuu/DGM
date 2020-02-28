clear all;
close all;
tic
a=0;
b=1;
resolution=50;
H=(b-a)/resolution;
t=a:H:b;

%Lets solve ODE du/dt=f(t,u)=3*t-u, u(0)=1
%Galerkin discretization
N=3;
P=4;%degree of approximation
h=(b-a)/N;
X=a:h:b;


    LPfast=cell(11,1);
    LPfast{1} = @(x) 1+x-x;
    LPfast{2} = @(x) x;
    LPfast{3} = @(x) 1/2*(3*x.^2-1);
    LPfast{4} = @(x) 1/2*(5*x.^3-3.*x);
    LPfast{5} = @(x) 1/8*(35*x.^4-30*x.^2+3);
    LPfast{6} = @(x) 1/8*(63*x.^5-70*x.^3+15*x);
    LPfast{7} = @(x) 1/16*(231*x.^6-315*x.^4+105*x.^2-5);
    LPfast{8} = @(x) 1/16*(429*x.^7-693*x.^5+315*x.^3-35*x);
    LPfast{9} = @(x) 1/128*(6435*x.^8-12012*x.^6+6930*x.^4-1260*x.^2+35);
    LPfast{10} = @(x) 1/128*(12155*x.^9-25740*x.^7+18018*x.^5-4620*x.^3+315*x);
    LPfast{11} = @(x) 1/256*(46189*x.^10-109395*x.^8+90090*x.^6-30030*x.^4+3465*x.^2-63);

LPderivative=cell(11,1);
LPder{1} = @(x) 0;
LPder{2} = @(x) 1;
LPder{3} = @(x) 3*x;
LPder{4} = @(x) 1/2*(15*x.^2-3);
LPder{5} = @(x) 1/8*(140*x.^3-60*x);
LPder{6} = @(x) 1/8*(315*x.^4-210*x.^2+15);
LPder{7} = @(x) 1/16*(1386*x.^5-1260*x.^3+210*x);
LPder{8} = @(x) 1/16*(3003*x.^6-3465*x.^4+945*x.^2-35);
%LPder{9} = @(x) 1/128*();
%LPder{10} = @(x)
% LPder{11} = @(x)



selfStencil=ones(P+1);
selfStencil=kron(eye(N),selfStencil);
upwindStencil=ones(P+1);
upwindStencil(1:2:P+1,:)=-1;
temporal = zeros(N);
temporal(2:N,1:N-1)=eye(N-1);
upwindStencil=kron(temporal, upwindStencil);
%Numerical flux matrix F, FU-initial condition is the flux
F=selfStencil + upwindStencil


legendreDiagonal=diag(2./(2*(0:P)+1));
deltaXs=diag(X(2:N+1)./2-X(1:N)./2);
D=kron(deltaXs, legendreDiagonal);

%For the calculation of mass matrix
quadWeights = [0.0812743883615744, 0.1806481606948574, 0.2606106964029354, 0.3123470770400029, 0.3302393550012598, ...
            0.3123470770400029, 0.2606106964029354, 0.1806481606948574, 0.0812743883615744];
        
quadPoints = [-0.9681602395076261, -0.8360311073266358, -0.6133714327005904, -0.3242534234038089, 0, ...
            0.3242534234038089, 0.6133714327005904, 0.8360311073266358, 0.9681602395076261];  


massStencil=zeros(P+1);
for p=0:P
    for q=0:P
        massStencil(q+1,p+1)=sum(quadWeights.*(LPfast{p+1}(quadPoints) .* LPder{q+1}(quadPoints) ));
        
    end
end    
M=kron(eye(N),massStencil);


Q=zeros((P+1)*N,1);
for e=1:N
    for q=0:P
        T=( X(e).*(1-quadPoints)./2 + X(e+1)*(1+quadPoints)./2 );
        
        Q((e-1)*(P+1)+q+1)=sum(...
            quadWeights.*(3.* T ...
            .*LPfast{q+1}(quadPoints) ).*(X(e+1)-X(e))./2 ...
            );
    end
end    
Q(1:P+1)=Q(1:P+1)-((-1).^(1:P+1))';
U=(F+D-M)\Q;

U=reshape(U,P+1,N);

Y=plotLP(t,X,U);

toc
hold on
plot(t,Y)


function y = f(t,u)
    y=3*t-u;
    
end    