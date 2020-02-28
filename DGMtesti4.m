%Uses Euler method
clear all; close all;
t=0:0.0005:5.4;
startT=t(1);
endT=t(end);
nT=length(t)
deltaT=(endT-startT)/(nT-1);

a=0;
b=10;

N=10;
P=6;
epsilon=0.00001;

x=a:0.001:b;
X=a:(b-a)/N:b;

iM=diag(kron(ones(1,N),(2*(0:P)+1)/2));
iJ=diag(kron(2./diff(X),ones(1,P+1)));
[currentUpper, nextLower, currentLower, previousUpper, coefficients] = centralFluxMatrices(N,P);

resultsU=zeros(N*(P+1),nT);
resultsA=zeros(N*(P+1),nT);
%IC=reshape(ApprLP(@initialA,X,P),N*(P+1),1);
IC=zeros(N*(P+1),1);
IC(1:P+1:end)=10.6;
%resultsA=kron(IC,ones(1,nT));
resultsA(:,1)=IC;
IC=zeros(N*(P+1),1);
IC(1:P+1:end)=1;
resultsU(:,1)=IC;
resultsA(:,1)=reshape(ApprLP(@initialA,X,P),N*(P+1),1);
%resultsU(:,1)=reshape(ApprLP(@initialU,X,P),N*(P+1),1);

%plot(x,plotLP(x,X,reshape(resultA(:,1),P+1,N) ));
% AL=1.2;
% AR=1;
[qW,qP]=GaussLegendreQuad(21);
LP=LegendreMatrix(qP);
LPminus=LegendreMatrix(qP-epsilon);
LPplus=LegendreMatrix(qP+epsilon);
LPder=LegendreDerivativeMatrix(qP);

pext=1000;
beta=2866;
a0=10.6;
%rho=0.001021;
rho=10.21;
MemoryA=zeros(N*(P+1),4);
MemoryU=zeros(N*(P+1),4);


tic
for i=1:nT
    RA=resultsA(:,i);
    RU=resultsU(:,i);
    
%     disp(num2str(i))
%     disp(currentUpper*RA)
%     disp(currentUpper*RU)
    
    resultsA(:,i+1)=RA + deltaT*RHS(N,RA,RU,@F1,@S1,qW,qP,LPminus,LPplus,LP,epsilon,bcA(i,deltaT),bcU(i,deltaT),10.6,0,beta,a0,rho, ...
        currentUpper, nextLower, currentLower, previousUpper,iM,iJ,LPder);
    resultsU(:,i+1)=RU + deltaT*RHS(N,RA,RU,@(a,u) F2(a,u,pext,beta,a0,rho),@(a,u) S2(a,u,rho),qW,qP,LPminus,LPplus,LP,epsilon,bcA(i,deltaT),bcU(i,deltaT),10.6,0,beta,a0,rho, ...
        currentUpper, nextLower, currentLower, previousUpper,iM,iJ,LPder);
    
    if i/floor(nT/10)==floor(i/floor(nT/10))
        fprintf('%s','*');
    end
end    
disp("Euler: " + toc)

% tic
% for i=1:nT
%     
% %Runge-kutta RK4
% %k1 = h*f(y(n))
% %k2 = h*f(y(n)+k1/2)
% %k3 = h*f(y(n)+k2/2)
% %k4 = h*f(y(n)+k3)
% %y(n+1)=y(n)+1/6*(k1+2*k2+2*k3+k4)
%     RA=resultsA(:,i);
%     RU=resultsU(:,i);
%     
%     %disp(num2str(i))
%     
%     k1a=deltaT*RHS(N,RA,RU,@F1,@S1,qW,qP,LPminus,LPplus,LP,epsilon,bcA(i,deltaT),bcU(i,deltaT),10.6,0,beta,a0,rho, ...
%         currentUpper, nextLower, currentLower, previousUpper,iM,iJ,LPder);
%     
%     k1u=deltaT*RHS(N,RA,RU,@(a,u) F2(a,u,pecurrentLowerxt,beta,a0,rho),@(a,u) S2(a,u,rho),qW,qP,LPminus,LPplus,LP,epsilon,bcA(i,deltaT),bcU(i,deltaT),10.6,0,beta,a0,rho, ...
%         currentUpper, nextLower, , previousUpper,iM,iJ,LPder);
%     
%     
%     k2a=deltaT*RHS(N,RA+k1a/2,RU+k1u/2,@F1,@S1,qW,qP,LPminus,LPplus,LP,epsilon,bcA(i,deltaT),bcU(i,deltaT),10.6,0,beta,a0,rho, ...
%         currentUpper, nextLower, currentLower, previousUpper,iM,iJ,LPder);
%     
%     k2u=deltaT*RHS(N,RA+k1a/2,RU+k1u/2,@(a,u) F2(a,u,pext,beta,a0,rho),@(a,u) S2(a,u,rho),qW,qP,LPminus,LPplus,LP,epsilon,bcA(i,deltaT),bcU(i,deltaT),10.6,0,beta,a0,rho, ...
%         currentUpper, nextLower, currentLower, previousUpper,iM,iJ,LPder);
%     
%     
%     k3a=deltaT*RHS(N,RA+k2a/2,RU+k2u/2,@F1,@S1,qW,qP,LPminus,LPplus,LP,epsilon,bcA(i,deltaT),bcU(i,deltaT),10.6,0,beta,a0,rho, ...
%         currentUpper, nextLower, currentLower, previousUpper,iM,iJ,LPder);
%     
%     k3u=deltaT*RHS(N,RA+k2a/2,RU+k2u/2,@(a,u) F2(a,u,pext,beta,a0,rho),@(a,u) S2(a,u,rho),qW,qP,LPminus,LPplus,LP,epsilon,bcA(i,deltaT),bcU(i,deltaT),10.6,0,beta,a0,rho, ...
%         currentUpper, nextLower, currentLower, previousUpper,iM,iJ,LPder);
%     
%     
%     k4a=deltaT*RHS(N,RA+k3a,RU+k3u,@F1,@S1,qW,qP,LPminus,LPplus,LP,epsilon,bcA(i,deltaT),bcU(i,deltaT),10.6,0,beta,a0,rho, ...
%         currentUpper, nextLower, currentLower, previousUpper,iM,iJ,LPder);
%     
%     k4u=deltaT*RHS(N,RA+k3a,RU+k3u,@(a,u) F2(a,u,pext,beta,a0,rho),@(a,u) S2(a,u,rho),qW,qP,LPminus,LPplus,LP,epsilon,bcA(i,deltaT),bcU(i,deltaT),10.6,0,beta,a0,rho, ...
%         currentUpper, nextLower, currentLower, previousUpper,iM,iJ,LPder);
%     
%     
%     resultsA(:,i+1)=RA+1/6*(k1a+2*k2a+2*k3a+k4a);
%     resultsU(:,i+1)=RU+1/6*(k1u+2*k2u+2*k3u+k4u);
% 
%     if i/floor(nT/10)==floor(i/floor(nT/10))
%         fprintf('%s','*');
%     end     
% end    
% disp("Runge-Kutta. " + toc)

fps=30;
frameT=floor(nT/(fps*(endT-startT)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%Plotting%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:frameT:nT
    tic
    A=reshape(resultsA(:,i),P+1,N);
    U=reshape(resultsU(:,i),P+1,N);
    
    Y=plotLP(x,X,A);
    Z=plotLP(x,X,U);
    plot(x,Y,'r');
    hold on
    plot(x,Z,'b');
    hold off
    axis([a b -0.5 15]);
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
    
    Y=plotLP(x,X,A);
    Z=plotLP(x,X,U);
    plot(x,Y,'r');
    hold on
    plot(x,Z,'b');
    axis([a b -0.5 15]);
    %text(1.02,b+0.1,'Time');
    %text(1.02,1,num2str((i-1)*(endT-startT)*1/nT));
    %disp(1/fps-toc)
    %pauseT=max( 1/fps-toc,0 )
end  


%disp(norm(result))


function y=F1(a,u)
    y=a.*u;
end

function y=F2(a,u,pext,beta,a0,rho)
    %a and u are the coefficient matrices of A and U
    %pext is the external pressure
    %beta is a elasticity parameter
    %a0 is the initial cross-sectional area
    %rho is the density of blood.
    y=u.*u/2+(pext+(beta./a0).*(sqrt(a)-sqrt(a0)))/rho;
    
    %Computational modelling of 1D blood...
    %y=u.*u/2+(pext+(beta)./a0.*(sqrt(a)-sqrt(a0)))/rho;
end  

function y=S1(a,u)
    y=zeros(size(a));
end

function y=S2(a,u,rho)
    y=-22*pi*0.0025*u ./(rho*a);
    %y=zeros(size(a));
end

function y=initialA(x)
    %y=5*sin(x)+20;
    %y=20;
    %y=0.106+0.05;
    %y=10.6+1*sin(x).^2;
    y=10.6;
end

function y=initialU(x)
    y=1;
end

% function y=W1(a,u,beta,a0,rho)
%     u0=0; 
%     gamma=sqrt(beta./(2*rho*a0));
%     y=u-u0+4*(gamma.*(nthroot(a,4) - nthroot(a0,4)));
% end
% 
% function y=W2(a,u,beta,a0,rho)
%     u0=0;
%     gamma=sqrt(beta./(2*rho*a0));
%     y=u-u0-4*(gamma.*(nthroot(a,4) - nthroot(a0,4)));
% end
% 
% function y=Aupwinded(a1,u1,a2,u2,beta,a0,rho)
%     y=nthroot( (sqrt(rho.*a0./(beta*32)).*(W1(a1,u1,beta,a0,rho)-W2(a2,u2,beta,a0,rho)) + nthroot(a0,4))  ,4);
% end
% 
% function y=Uupwinded(a1,u1,a2,u2,beta,a0,rho)
%     y=0.5*(W1(a1,u1,beta,a0,rho)+W2(a2,u2,beta,a0,rho));
% 
% end

function y=bcA(i,deltaT)
    %y=sin(10*deltaT*(i-1));
    %y=5*sin(6*deltaT*(i-1))+20;
    %y=20+25*exp(-0.5*((deltaT*(i-1)-1)/0.2).^2);
    %y=0.106+0.1*sin((deltaT*(i-1))).^2;
    %y=10.6;
    y=10.6+1.5+1.5*cos(4*deltaT*(i-1));
    %y=floor(i/100)-2*floor(i/200);
    %y=10+i*deltaT;
    %y=0;
end 
function y=bcU(i,deltaT)
    %y=0.5*sin(2*deltaT*(i-1))+1;
    %y=15*sin(10*deltaT*(i-1))+15;
    %y=floor(i/100)-2*floor(i/200);
    %y=1+deltaT*i;
    y=1;
    %y=1.5-0.5*cos(2*deltaT*(i-1));
end 

