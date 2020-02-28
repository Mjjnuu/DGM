clear all; close all;

betaL=304;
betaR=304;
AL=0.116;
AR=0.106;
UL=1.1;
UR=1.2;
rho=1050;
AL0=0.106;
AR0=0.106;

iterations=10;

time1=0;
for i=1:100
    tic
    y=NewtonRhapson(@(x) F(x,AL,UL,AR,UR,AL0,AR0,betaL,betaR,rho),4,[AL UL AR UR],4,iterations);
    time1=time1+toc;
end
%fprintf('%6.12f\n',y)
disp("   Normal: " + time1)

time2=0;
for i=1:100
    tic
    y=NewtonRhapson(@(x) fastF(x,AL,UL,AR,UR,AL0,AR0,betaL,betaR,rho),4,[AL UL AR UR],4,iterations);
    time2=time2+toc;
end
%fprintf('%6.12f\n',y)
disp("   Fast: " + time2)

time3=0;
for i=1:100
    tic
    y=FastNewtonRhapson(@(x) fastF(x,AL,UL,AR,UR,AL0,AR0,betaL,betaR,rho), [AL UL AR UR], iterations, AL,UL,AR,UR,AL0,AR0,betaL,betaR,rho);
    time3=time3+toc;
end
%fprintf('%6.12f\n',y)
disp("   Fastest: " + time3)

disp("   Ratio. " + time1/time3)
disp(" ")
FastNewtonRhapson(@(x) fastF(x,AL,UL,AR,UR,AL0,AR0,betaL,betaR,rho), [AL UL AR UR], iterations, AL,UL,AR,UR,AL0,AR0,betaL,betaR,rho)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function data=NewtonRhapson(f,m,x0,n,iterations)
%     data=[0,x0,0];
%     for i=1:iterations
%         x1=x0'-numJacob(f,m,x0,n)\f(x0);
%         error=norm(x1'-x0);
%         x0=x1';
%         data=[data;[i x0 error]];
%     end    
% end

function y=NewtonRhapson(f,m,x0,n,iterations)
    for i=1:iterations
        x1=x0'-numJacob(f,m,x0,n)\f(x0);
        x0=x1';
    end 
    y=x0;
end

function y=FastNewtonRhapson(f,x0,iterations,AL,UL,AR,UR,AL0,AR0,betaL,betaR,rho)
    for i=1:iterations
        x1=x0'-FastNumJacob(x0,AL,UL,AR,UR,AL0,AR0,betaL,betaR,rho)\f(x0);
        x0=x1';
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

function Jf=FastNumJacob(x,AL,UL,AR,UR,AL0,AR0,betaL,betaR,rho)
    Jf=zeros(4);
    Jf(1,1)=-sqrt(betaL/(2*rho*AL0))*x(1).^(-0.75);
    Jf(1,2)=-1;
    Jf(2,3)=sqrt(betaR/(2*rho*AR0))*x(3).^(-0.75);
    Jf(2,4)=-1;
    Jf(3,1)=x(2);
    Jf(3,2)=x(1);
    Jf(3,3)=-x(4);
    Jf(3,4)=-x(3);
    Jf(4,1)=0.5*betaL/AL0*x(1).^(-0.5);
    Jf(4,2)=rho*x(2);
    Jf(4,3)=-0.5*betaR/AR0*x(3).^(-0.5);
    Jf(4,4)=-rho*x(4);
    
end

%variables x(1)=ALu, x(2)=ULu, x(3)=ARu and x(4)=URu
function y=F(x,AL,UL,AR,UR,AL0,AR0,betaL,betaR,rho)
    y(1,1)=UL-x(2)+4*sqrt(betaL/(2*rho.*AL0))*(nthroot(AL,4)-nthroot(x(1),4));
    y(2,1)=UR-x(4)-4*sqrt(betaR/(2*rho.*AR0))*(nthroot(AR,4)-nthroot(x(3),4));
    y(3,1)=x(1).*x(2)-x(3).*x(4);
    y(4,1)=rho*0.5*(x(2).^2-x(4).^2)+betaL/AL0*(sqrt(x(1))-sqrt(AL0))-betaR/AR0*(sqrt(x(3))-sqrt(AR0));
end

function y=fastF(x,AL,UL,AR,UR,AL0,AR0,betaL,betaR,rho)
    y(1,1)=UL-x(2)+4*sqrt(betaL/(2*rho.*AL0))*(AL.^0.25-x(1).^0.25);
    y(2,1)=UR-x(4)-4*sqrt(betaR/(2*rho.*AR0))*(AR.^0.25-x(3).^0.25);
    y(3,1)=x(1).*x(2)-x(3).*x(4);
    y(4,1)=rho*0.5*(x(2).^2-x(4).^2)+betaL/AL0*(sqrt(x(1))-sqrt(AL0))-betaR/AR0*(sqrt(x(3))-sqrt(AR0));
end


function y=accuracyTest(x,AL,UL,AR,UR,AL0,AR0,betaL,betaR,rho)
    y(1,1)=UL+4*sqrt(betaL/(2*rho*AL0))*(nthroot(AL,4)-nthroot(AL0,4));
    y(1,2)=x(2)+4*sqrt(betaL/(2*rho*AL0))*(nthroot(x(1),4)-nthroot(AL0,4));
    y(2,1)=UR-4*sqrt(betaR/(2*rho*AR0))*(nthroot(AR,4)-nthroot(AR0,4));
    y(2,2)=x(4)-4*sqrt(betaR/(2*rho*AR0))*(nthroot(x(3),4)-nthroot(AR0,4));
    y(3,1)=x(1).*x(2);
    y(3,2)=x(3).*x(4);
    y(4,1)=rho*(x(2).^2)/2+betaL/AL0*(sqrt(x(1)) - sqrt(AL0));
    y(4,2)=rho*(x(4).^2)/2+betaR/AR0*(sqrt(x(3)) - sqrt(AR0));
end






