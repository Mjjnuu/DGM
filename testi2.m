clear all; close all;
%F=@(x) -x+2*x.^2+2*x.^3;
%dF=@(x) -1+4*x+6*x.^2;
F=@(t) t.^2;
dF=@(t) 2*t;
id=@(x) x;
u=@(t) exp(t-t.^2);
du=@(t) (1-2*t).*exp(t-t.^2);

P=7;
N=1;
x=-1:0.001:1;
X=-1:2/N:1;
U=ApprLP(u,X,P);
U=reshape(U,N*(P+1),1);
size(U);
[qW,qP]=GaussLegendreQuad(21);

epsilon=0.0000001;

LP=LegendreMatrix(qP);
LPminus=LegendreMatrix(qP-epsilon);
LPplus=LegendreMatrix(qP+epsilon);

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

plot(x,dF(u(x)).*du(x),'k');
hold on
plot(x,F(u(x)),'r')
plot(x,F(plotLP(x,X,(reshape(U,P+1,N)))))
scatter(qP,F(reshape(U,P+1,N)'*LP(1:P+1,:)))

symmetricDifference=F(reshape(U,P+1,N)'*LPplus(1:P+1,:))-F(reshape(U,P+1,N)'*LPminus(1:P+1,:));
%disp(size(symmetricDifference));
symmetricDifference=symmetricDifference/(2*epsilon);

scatter(qP,symmetricDifference)


testaus=zeros(N*(P+1),1);

for e=1:N
    for q=0:P
        product=@(t) 2*exp(t-t.^2).*(1-2*t).*exp(t-t.^2).*LPfast{q+1}((2*t-(X(e)+X(e+1)))/(X(e+1)-X(e)));
        testaus((e-1)*(P+1)+q+1)=integral(product,X(e),X(e+1));
    end    
end    

testaus=testaus';
approximation=derivativeIntegral(N,U,F,qW,qP,LPminus,LPplus,LP,epsilon)';
size(approximation);
size(testaus);
error=approximation-testaus;
L2error=sqrt(sum(error.*error));
L2testaus=sqrt(sum(testaus.*testaus));
L2approximation=sqrt(sum(approximation.*approximation));
disp("L2 norm of the error: " + L2error);
disp("L2 norm of the testaus: "+L2testaus)
disp("L2 norm of the approximation: "+L2approximation)
disp("max error: " + max(abs(error)));

disp("Relative error: " + L2error/L2testaus)




