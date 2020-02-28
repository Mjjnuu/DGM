A1=0.55;
A2=0.126;
A3=0.158;
U1=6.5;
U2=7.4;
U3=9;
A10=0.51;
A20=0.106;
A30=0.145;
beta1=466;
beta2=2866;
beta3=2246;
rho=1021;

fprintf('%2.16f\n',1/sqrt(5.^3))
fprintf('%2.16f\n',5.^(-1.5))

tic
parfor i=1:1000
    data=NewtonRhapson(@(x) splitJunction(x,A1,U1,A2,U2,A3,U3,A10,A20,A30,beta1,beta2,beta3,rho),6,[A1 U1 A2 U2 A3 U3],6,6);
end    
disp("parfor: " + toc)

tic
for i=1:1000
    data=NewtonRhapson(@(x) splitJunction(x,A1,U1,A2,U2,A3,U3,A10,A20,A30,beta1,beta2,beta3,rho),6,[A1 U1 A2 U2 A3 U3],6,6);
end    
disp("for: " + toc)


time=0;
for i=1:1000
   tic
   1/sqrt(5.^3);
   time=time+toc;
end    
disp("sqrt: " + time)

time=0;
for i=1:1000
   tic
   5.^(-1.5);
   time=time+toc;
end    
disp("power: " + time)



fprintf('%2.16f  %2.16f\n', nthroot(1000,4), 1000.^(1/4))

% tic
% data=NewtonRhapson(@(x) splitJunction(x,A1,U1,A2,U2,A3,U3,A10,A20,A30,beta1,beta2,beta3,rho),6,[A1 U1 A2 U2 A3 U3],6,8)
% toc
time=0;
for i=1:100
    tic
    data=NewtonRhapson(@(x) splitJunction(x,A1,U1,A2,U2,A3,U3,A10,A20,A30,beta1,beta2,beta3,rho),6,[A1 U1 A2 U2 A3 U3],6,6);
    time=time+toc;
end 
disp(time/100)
% tic
% x=fsolve(@(x) splitJunction(x,A1,U1,A2,U2,A3,U3,A10,A20,A30,beta1,beta2,beta3,rho),[A1 U1 A2 U2 A3 U3])
% toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y=NewtonRhapson(f,m,x0,n,iterations)
    %data=[0,x0,0];
    for i=1:iterations
        x1=x0'-numJacob(f,m,x0,n)\f(x0);
        %error=norm(x1'-x0);
        x0=x1';
        %data=[data;[i x0 error]];
    end 
    y=x0;
end


function Jf=numJacob(f,m,x,n)
    Jf=ones(m,n);
    h=0.0000001;
    for j=1:n
        e=zeros(1,n); e(j)=1;
        Jf(:,j)=( f(x+e*h)-f(x-e*h) )/(2*h);
        
    end    
end


function y=splitJunction(x,A1,U1,A2,U2,A3,U3,A10,A20,A30,beta1,beta2,beta3,rho)
    y(1,1)=U1-x(2)+4*sqrt(beta1/(2*rho.*A10))*(nthroot(A1,4)-nthroot(x(1),4));
    y(2,1)=U2-x(4)-4*sqrt(beta2/(2*rho.*A20))*(nthroot(A2,4)-nthroot(x(3),4));
    y(3,1)=U3-x(6)-4*sqrt(beta3/(2*rho.*A30))*(nthroot(A3,4)-nthroot(x(5),4));
    
    y(4,1)=x(1).*x(2)-x(3).*x(4)-x(5).*x(6);
    y(5,1)=rho*0.5*(x(2).^2-x(4).^2)+beta1/A10*(sqrt(x(1))-sqrt(A10))-beta2/A20*(sqrt(x(3))-sqrt(A20));
    y(6,1)=rho*0.5*(x(2).^2-x(6).^2)+beta1/A10*(sqrt(x(1))-sqrt(A10))-beta3/A30*(sqrt(x(5))-sqrt(A30));

end

