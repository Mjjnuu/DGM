clear all; close all;
t=0:0.0001:5;
startT=t(1);
endT=t(end);
nT=length(t)
deltaT=(endT-startT)/(nT-1);

a=0;
b=10;
aa=10;
bb=20;
aaa=20;
bbb=25;

N=10;
P=3;
epsilon=0.000001;

pext=77310; %773.1Pa=kg*m^-1*s^-2=100*kg*cm^-1*s^-2
%pext=1.021*10^(-3);
%pext=1021;

betaL=303; betaR=330; beta3=350;
AL0=0.106; AR0=0.09; A30=0.08;
%rho=0.001021;
%rho=10.21;
rho=1;
%rho=1.021*10^(-3);
%rho=1.035*10^(-3);

x=a:0.005:b;
X=a:(b-a)/N:b;
xx=aa:0.005:bb;
XX=aa:(bb-aa)/N:bb;
xxx=aaa:0.005:bbb;
XXX=aaa:(bbb-aaa)/N:bbb;

iM=diag(kron(ones(1,N),(2*(0:P)+1)/2));
iJ=diag(kron(2./diff(X),ones(1,P+1)));
iJJ=diag(kron(2./diff(XX),ones(1,P+1)));
iJJJ=diag(kron(2./diff(XXX),ones(1,P+1)));
[currentUpper, nextLower, currentLower, previousUpper, coefficients] = centralFluxMatrices(N,P);

resultsU=zeros(N*(P+1),nT);   resultsUU=zeros(N*(P+1),nT);   resultsUUU=zeros(N*(P+1),nT);
resultsA=zeros(N*(P+1),nT);   resultsAA=zeros(N*(P+1),nT);   resultsAAA=zeros(N*(P+1),nT);
%IC=reshape(ApprLP(@initialA,X,P),N*(P+1),1);
IC=zeros(N*(P+1),1);
IC(1:P+1:end)=AL0;  resultsA(:,1)=IC;
IC(1:P+1:end)=AR0;  resultsAA(:,1)=IC;
IC(1:P+1:end)=A30;  resultsAAA(:,1)=IC;
IC=zeros(N*(P+1),1);
IC(1:P+1:end)=0;
resultsU(:,1)=IC;  resultsUU(:,1)=IC; resultsUUU(:,1)=IC;%Both flow rates are initially zero
%resultsA(:,1)=reshape(ApprLP(@initialA,X,P),N*(P+1),1);
%resultsU(:,1)=reshape(ApprLP(@initialU,X,P),N*(P+1),1);

[qW,qP]=GaussLegendreQuad(9);
LP=LegendreMatrix(qP);
LPminus=LegendreMatrix(qP-epsilon);
LPplus=LegendreMatrix(qP+epsilon);
LPder=LegendreDerivativeMatrix(qP);

MemoryA=zeros(N*(P+1),4);  MemoryAA=zeros(N*(P+1),4);   MemoryAAA=zeros(N*(P+1),4);
MemoryU=zeros(N*(P+1),4);  MemoryUU=zeros(N*(P+1),4);   MemoryUUU=zeros(N*(P+1),4);


disp("         |         |         |         |         |         |         |         |         |         |")
tic
for i=1:nT
    RA=resultsA(:,i);  RAA=resultsAA(:,i);  RAAA=resultsAAA(:,i);
    RU=resultsU(:,i);  RUU=resultsUU(:,i);  RUUU=resultsUUU(:,i);
    AL=ones(1,P+1)*RA((N-1)*(P+1)+1:N*(P+1));
    UL=ones(1,P+1)*RU((N-1)*(P+1)+1:N*(P+1));
    AR=(-1).^(0:P)*RAA((N-1)*(P+1)+1:N*(P+1));
    UR=(-1).^(0:P)*RUU((N-1)*(P+1)+1:N*(P+1));
    upwinded=FastNewtonRhapson([AL UL AR UR],10,AL,UL,AR,UR,AL0,AR0,betaL,betaR,rho);
    
    AL=ones(1,P+1)*RAA((N-1)*(P+1)+1:N*(P+1));
    UL=ones(1,P+1)*RUU((N-1)*(P+1)+1:N*(P+1));
    AR=(-1).^(0:P)*RAAA((N-1)*(P+1)+1:N*(P+1));
    UR=(-1).^(0:P)*RUUU((N-1)*(P+1)+1:N*(P+1));
    UW=FastNewtonRhapson([AL UL AR UR],10,AL,UL,AR,UR,AR0,A30,betaR,beta3,rho);
    
    resultsA(:,i+1)=RA + deltaT*RHS(N,RA,RU,@F1,@S1,qW,qP,LPminus,LPplus,LP,epsilon,bcA(i,deltaT),bcU(i,deltaT),upwinded(1),upwinded(2),betaL,AL0,rho, ...
        currentUpper, nextLower, currentLower, previousUpper,iM,iJ,LPder);
    resultsU(:,i+1)=RU + deltaT*RHS(N,RA,RU,@(a,u) F2(a,u,pext,betaL,AL0,rho),@(a,u) S2(a,u,rho),qW,qP,LPminus,LPplus,LP,epsilon,bcA(i,deltaT),bcU(i,deltaT),upwinded(1),upwinded(2),betaL,AL0,rho, ...
        currentUpper, nextLower, currentLower, previousUpper,iM,iJ,LPder);
    
    
    resultsAA(:,i+1)=RAA + deltaT*RHS(N,RAA,RUU,@F1,@S1,qW,qP,LPminus,LPplus,LP,epsilon,upwinded(3),upwinded(4),UW(1),UW(2),betaR,AR0,rho, ...
        currentUpper, nextLower, currentLower, previousUpper,iM,iJJ,LPder);
    resultsUU(:,i+1)=RUU + deltaT*RHS(N,RAA,RUU,@(a,u) F2(a,u,pext,betaR,AR0,rho),...
        @(a,u) S2(a,u,rho),qW,qP,LPminus,LPplus,LP,epsilon,upwinded(3),upwinded(4),UW(1),UW(2),betaR,AR0,rho, ...
        currentUpper, nextLower, currentLower, previousUpper,iM,iJJ,LPder);
    
    
    resultsAAA(:,i+1)=RAAA + deltaT*RHS(N,RAAA,RUUU,@F1,@S1,qW,qP,LPminus,LPplus,LP,epsilon,UW(3),UW(4),A30,0,beta3,A30,rho, ...
        currentUpper, nextLower, currentLower, previousUpper,iM,iJJJ,LPder);
    resultsUUU(:,i+1)=RUUU + deltaT*RHS(N,RAAA,RUUU,@(a,u) F2(a,u,pext,beta3,A30,rho),...
        @(a,u) S2(a,u,rho),qW,qP,LPminus,LPplus,LP,epsilon,UW(3),UW(4),A30,0,beta3,A30,rho, ...
        currentUpper, nextLower, currentLower, previousUpper,iM,iJJJ,LPder);
    
    if i/floor(nT/100)==floor(i/floor(nT/100))
        fprintf('%s','*');
    end
end    
fprintf('%s\n', '');
disp("Euler: " + toc)


fps=24;
frameT=floor(nT/(fps*(endT-startT)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%Plotting%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
multiplier=1;
for i=1:frameT:nT
    tic
    A=reshape(resultsA(:,i),P+1,N);
    U=reshape(resultsU(:,i),P+1,N);
    AA=reshape(resultsAA(:,i),P+1,N);
    UU=reshape(resultsUU(:,i),P+1,N);
    AAA=reshape(resultsAAA(:,i),P+1,N);
    UUU=reshape(resultsUUU(:,i),P+1,N);
    
    Y=plotLP(x,X,A);
    Z=plotLP(x,X,U);
    YY=plotLP(xx,XX,AA);
    ZZ=plotLP(xx,XX,UU);
    YYY=plotLP(xxx,XXX,AAA);
    ZZZ=plotLP(xxx,XXX,UUU);
%     Y=multiplier*(betaL/AL0*(sqrt(Y)-sqrt(AL0))+pext);
%     YY=multiplier*(betaR/AR0*(sqrt(YY)-sqrt(AR0))+pext);
%     YYY=multiplier*(beta3/A30*(sqrt(YYY)-sqrt(A30))+pext);
    
    plot([x xx xxx],[Y YY YYY],'r');
    hold on
    plot([x xx xxx],[Z ZZ ZZZ],'b');
    
    hold off
    axis([a bbb -1 20]);
    text(1.02,b+0.1,'Time');
    text(1.02,1,num2str((i-1)*(endT-startT)*1/nT));
    %disp(1/fps-toc)
    %pauseT=max( 1/fps-toc,0 )
    pause(max( 1/fps-toc,0 ));
end  

figure
for i=1:floor(nT/10):nT
    tic
    A=reshape(resultsA(:,i),P+1,N);
    U=reshape(resultsU(:,i),P+1,N);
    AA=reshape(resultsAA(:,i),P+1,N);
    UU=reshape(resultsUU(:,i),P+1,N);
    AAA=reshape(resultsAAA(:,i),P+1,N);
    UUU=reshape(resultsUUU(:,i),P+1,N);
    
    Y=plotLP(x,X,A);
    Z=plotLP(x,X,U);
    YY=plotLP(xx,XX,AA);
    ZZ=plotLP(xx,XX,UU);
    YYY=plotLP(xxx,XXX,AAA);
    ZZZ=plotLP(xxx,XXX,UUU);
%     Y=multiplier*(betaL/AL0*(sqrt(Y)-sqrt(AL0))+pext);
%     YY=multiplier*(betaR/AR0*(sqrt(YY)-sqrt(AR0))+pext);
%     YYY=multiplier*(beta3/A30*(sqrt(YYY)-sqrt(A30))+pext);

    
    plot([x xx xxx],[Y YY YYY],'r');
    hold on
    plot([x xx xxx],[Z ZZ ZZZ],'b');
    axis([a bbb -1 20]);
end  

%disp(norm(result))


function y=F1(a,u)
    y=a.*u;
end

function y=F2(a,u,pext,beta,a0,rho)
    y=u.*u/2+(pext+(beta./a0).*(sqrt(a)-sqrt(a0)))/rho;
    
    %Computational modelling of 1D blood...
    %y=u.*u/2+(pext+(beta)./a0.*(sqrt(a)-sqrt(a0)))/rho;
end  

function y=S1(a,u)
    y=zeros(size(a));
end

function y=S2(a,u,rho)
    y=-22*pi*0.0025*u ./(rho*a);
    %y=-22*pi*u ./(rho*a);
    %y=zeros(size(a));
end

function y=initialA(x)
    %y=10.6+1*sin(x).^2;
    y=0.106;
end

function y=initialU(x)
    y=0;
end

function y=bcA(i,deltaT)
    y=0.106;
    %y=10.6+1.5+1.5*sin(4*deltaT*(i-1));
    %y=floor(i/100)-2*floor(i/200);
    %y=10+i*deltaT;
    %y=0;
end 
function y=bcU(i,deltaT)
    %y=1+deltaT*i;
    %y=1;
    %y=1.5-0.5*cos(2*deltaT*(i-1));
    y=15*(1-cos(8*deltaT*(i-1)));
end 



function y=FastNewtonRhapson(x0,iterations,AL,UL,AR,UR,AL0,AR0,betaL,betaR,rho)
    for i=1:iterations
        x1=x0'-FastNumJacob(x0,AL0,AR0,betaL,betaR,rho)\fastF(x0,AL,UL,AR,UR,AL0,AR0,betaL,betaR,rho);
        x0=x1';
    end 
    y=x0;
end

function Jf=FastNumJacob(x,AL0,AR0,betaL,betaR,rho)
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

function y=fastF(x,AL,UL,AR,UR,AL0,AR0,betaL,betaR,rho)
    y(1,1)=UL-x(2)+4*sqrt(betaL/(2*rho.*AL0))*(AL.^0.25-x(1).^0.25);
    y(2,1)=UR-x(4)-4*sqrt(betaR/(2*rho.*AR0))*(AR.^0.25-x(3).^0.25);
    y(3,1)=x(1).*x(2)-x(3).*x(4);
    y(4,1)=rho*0.5*(x(2).^2-x(4).^2)+betaL/AL0*(sqrt(x(1))-sqrt(AL0))-betaR/AR0*(sqrt(x(3))-sqrt(AR0));
end