%Testing of the Gauss-Lobatto quadrature
%\int_{-1}^{1}f(x)dx = 2/n*(n-1)*(f(1)-f(-1)) + sum_{1=2}^{n-1} wi*f(xi)

%n=3, points xi = 0, +-1: wi = 4/3, 1/3
%n=4, points xi = +-sqrt(1/5),+-1; wi = 5/6, 1/6
%n=5, points xi = 0, +-sqrt(3/7), +-1; wi = 32/45, 49/90, 1/10
%n=6


%Gauss-Legendre
%n=6 points [-0.9324695142031521, -0.6612093864662645, -0.2386191860831969, 0.2386191860831969, 0.6612093864662645, 0.9324695142031521]
%    weights [0.1713244923791704, 0.3607615730481386, 0.4679139345726910, 0.4679139345726910, 0.3607615730481386, 0.1713244923791704]

%n=9 points [-0.9681602395076261, -0.8360311073266358, -0.6133714327005904, -0.3242534234038089, 0, ...
%   0.3242534234038089, 0.6133714327005904, 0.8360311073266358, 0.9681602395076261]
%
%    weghts [0.0812743883615744, 0.1806481606948574, 0.2606106964029354, 0.3123470770400029, 0.3302393550012598, ...
%             0.3123470770400029, 0.2606106964029354, 0.1806481606948574, 0.1806481606948574, 0.0812743883615744]
t=-1:0.0001:1;

f = @(x) exp(x).*sin(x)-cos(x)+sqrt(x.^2+sqrt(x.^2+sqrt(x.^2+1))) + exp(x).^sin(sin(cos(cos(x)))) - sqrt(x);
id = @(x) 1;
acc=exp(1)-exp(-1);
plot(t,f(t))

inttime=0;
for i=1:1000
    tic
    int=integral(f,-1,1,'RelTol',0.001,'AbsTol',1e-3);
    inttime = inttime + toc;
end  


gltime=0;
for i=1:1000
    tic
%     quadweights=[1/10, 49/90, 32/45, 49/90, 1/10];
%     quadpoints=[-1, -sqrt(3/7), 0, sqrt(3/7), 1];
%     fi=f(quadpoints);
%     gausslobattointegral=sum(fi.*quadweights);
    gausslobattointegral = FastL2norm(f,id);
    gltime= gltime + toc;
end    
  

disp(" ");
fprintf('%20s    %12.12f   %e\n', 'Gauss-Lobatto:' , gltime , gausslobattointegral-acc);
fprintf('%20s    %12.12f   %e\n', 'Integral:' , inttime , int-acc);
fprintf('%20s    %12.12f\n', 'Ratio:', inttime/gltime);
disp(" ");
