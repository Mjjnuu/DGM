clear all
n=100;
h=2/n;
t=-2:h:2;

%Approksimoidaan funktiota e^t-1.2sin(t)
f = @(x) exp(x)-1.9*sin(3*x)+1/15*x.^5;
%f = @(x) x.^5+2*x.^4 - x.^3+0.5*x.^2+x-1;

l0 = @(x) legendreP(0,x);
l1 = @(x) legendreP(1,x);
l2 = @(x) legendreP(2,x);
l3 = @(x) legendreP(3,x);
l4 = @(x) legendreP(4,x);
l5 = @(x) legendreP(5,x);
l6 = @(x) legendreP(6,x);
l7 = @(x) legendreP(7,x);
l8 = @(x) legendreP(8,x);
l9 = @(x) legendreP(9,x);

tic
F=zeros(5,1);%miksi t‰ss‰ oli (2,1)????
% F(1)=L2norm(l0,f,-1,1);
% F(2)=L2norm(l1,f,-1,1);
% F(3)=L2norm(l2,f,-1,1);
% F(4)=L2norm(l3,f,-1,1);
% F(5)=L2norm(l4,f,-1,1);

for i=1:10
    F(i)=LegendreProduct(f,i-1,-1,1);
end

d=2./(2*(0:9)+1);
M=diag(2./(2*(0:9)+1));
u=M\F;
A=[l0(t);l1(t);l2(t);l3(t);l4(t);l5(t);l6(t);l7(t);l8(t);l9(t)];
B=u'*A;  
%Tavallisia legendren polynomeja k‰ytet‰‰n t‰ss‰ plottaamiseen, joten 
%integrointirajat legendren tulossa ei kasvata kuvaajan tarkkuutta

disp("Time: " + toc)

plot (t,f(t));
hold on 
%plot (t,Y);
plot (t,B)
%L2-norm of the error:
disp("L2 error: " + sum(((f(t)-B).*h).^2))

%plot(t,LegendreTransform(7,t))


% function norm = L2norm(a,b)
%     product = @(x) a(x).*b(x);
%     norm = integral(product,-1,1);
% end    


