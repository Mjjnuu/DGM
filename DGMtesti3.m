close all 
clear all

%9th order Gauss-Legendre
% qW = [0.0812743883615744; 0.1806481606948574; 0.2606106964029354; ...
%     0.3123470770400029; 0.3302393550012598; 0.3123470770400029; ...
%     0.2606106964029354; 0.1806481606948574; 0.0812743883615744];
% 
% qP = [-0.9681602395076261, -0.8360311073266358, -0.6133714327005904, ...
%     -0.3242534234038089, 0,0.3242534234038089, ...
%     0.6133714327005904, 0.8360311073266358, 0.9681602395076261];

%11th order Gauss-Legendre
% qW=[0.0556685671161737; 0.1255803694649046; 0.1862902109277343; ...
%     0.2331937645919905; 0.2628045445102467; 0.2729250867779006; ...
%     0.2628045445102467; 0.2331937645919905; 0.1862902109277343; ...
%     0.1255803694649046; 0.0556685671161737;];
% 
% qP=[-0.9782286581460570, -0.8870625997680953, -0.7301520055740494 ...
%     -0.5190961292068118, -0.2695431559523450, 0.0000000000000000 ...
%     0.2695431559523450, 0.5190961292068118, 0.7301520055740494 ...
%     0.8870625997680953, 0.9782286581460570];

%21th order Gauss-Legendre
qW=[0.0160172282577743; 0.0369537897708525; 0.0571344254268572; ...
    0.0761001136283793; 0.0934444234560339; 0.1087972991671484; ...
    0.1218314160537285; 0.1322689386333375; 0.1398873947910731; ...
    0.1445244039899700; 0.1460811336496904; 0.1445244039899700; ...
    0.1398873947910731; 0.1322689386333375; 0.1218314160537285; ...
    0.1087972991671484; 0.0934444234560339; 0.0761001136283793; ...
    0.0571344254268572; 0.0369537897708525; 0.0160172282577743];

qP=[-0.9937521706203895, -0.9672268385663063, -0.9200993341504008, ...
    -0.8533633645833173, -0.7684399634756779, -0.6671388041974123, ...
    -0.5516188358872198, -0.4243421202074388, -0.2880213168024011, ...
    -0.1455618541608951, 0.0000000000000000, 0.1455618541608951, ...
    0.2880213168024011, 0.4243421202074388, 0.5516188358872198, ...
    0.6671388041974123, 0.7684399634756779, 0.8533633645833173, ...
    0.9200993341504008, 0.9672268385663063, 0.9937521706203895];
    

    %4th order Gauss-Legendre
%     qW=[0.3478548451374538; 0.6521451548625461; 0.6521451548625461; 0.3478548451374538];
%     qP=[-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526];
% LP=zeros(11,length(qW));
% LP(1,:)=1;
% LP(2,:)=qP;
% LP(3,:)=1/2*(3*qP.^2-1);
% LP(4,:)=1/2*(5*qP.^3-3.*qP);
% LP(5,:)=1/8*(35*qP.^4-30*qP.^2+3);
% LP(6,:)=1/8*(63*qP.^5-70*qP.^3+15*qP);
% LP(7,:)=1/16*(231*qP.^6-315*qP.^4 ...
%     +105*qP.^2-5);
% 
% LP(8,:)=1/16*(429*qP.^7-693*qP.^5 ...
%     +315*qP.^3-35*qP);
% 
% LP(9,:)=1/128*(6435*qP.^8-12012*qP.^6 ...
%     +6930*qP.^4-1260*qP.^2+35);
% 
% LP(10,:)=1/128*(12155*qP.^9-25740*qP.^7 ...
%     +18018*qP.^5-4620*qP.^3+315*qP);
% 
% LP(11,:)=1/256*(46189*qP.^10-109395*qP.^8 ...
%     +90090*qP.^6-30030*qP.^4+3465*qP.^2-63);

LP=LegendreMatrix(qP);
LPder=LegendreDerivativeMatrix(qP);

% LPder=zeros(11,length(qW));    
% LPder(1,:) = 0;
% LPder(2,:) = 1;
% LPder(3,:) = 3*qP;
% LPder(4,:) = 1/2*(15*qP.^2-3);
% LPder(5,:) = 1/8*(140*qP.^3-60*qP);
% LPder(6,:) = 1/8*(315*qP.^4-210*qP.^2+15);
% LPder(7,:) = 21/8*(33*qP.^5-30*qP.^3+5*qP);
% LPder(8,:) = 7/16*(429*qP.^6-495*qP.^4+135*qP.^2-5);
% LPder(9,:) = 9/16*(715*qP.^7-1001*qP.^5+385*qP.^3-35*qP);
% LPder(10,:) = 45/128*(2431*qP.^8-4004*qP.^6+2002*qP.^4-308*qP.^2+7);
% LPder(11,:) = 55/128*(4199*qP.^9-7956*qP.^7+4914*qP.^5-1092*qP.^3+63);

    
% LPder=zeros(11,length(qW));    
% LPder(1,:) = 0;
% LPder(2,:) = 1;
% LPder(3,:) = 3*qP;
% LPder(4,:) = 1/2*(15*qP.^2-3);
% LPder(5,:) = 1/8*(140*qP.^3-60*qP);
% LPder(6,:) = 1/8*(315*qP.^4-210*qP.^2+15);
% LPder(7,:) = 1/16*(1386*qP.^5-1260*qP.^3+210*qP);
% LPder(8,:) = 1/16*(3003*qP.^6-3465*qP.^4+945*qP.^2-35);
% LPder(9,:) = 1/128*(51480*qP.^7-72072*qP.^5+27720*qP.^3-3780*qP);
% LPder(10,:) = 1/128*(109395*qP.^8-180180*qP.^6+90090*qP.^4-13860*qP.^2+315);
% LPder(11,:) = 1/256*(461890*qP.^9-875160*qP.^7+540540*qP.^5-120120*qP.^3+6930*qP);

%The goal is to solve a hyperbolic partial differential equation of the 
%form du/dt+df(u)/dx=s(u), by using the discontinuous Galerkin method and 
%explicit 4th order Runge-Kutta time integration.
%
%The spatial domain [a,b] is partitioned to N elements of equal length, and
%the polynomial approximation space contains polynomials of degree P or
%lower.
a=0;
b=2;
resolution=1000;
h=(b-a)/resolution;
x=a:h:b;

N=10;
P=5;

deltaX=(b-a)/N;
X=a:deltaX:b;

Mass=diag(kron(ones(1,N),2./(2*(0:P)+1)));
Jacobian=diag(kron(diff(X)./2,ones(1,P+1)));

iM=diag(kron(ones(1,N),(2*(0:P)+1)/2));
iJ=diag(kron(2./diff(X),ones(1,P+1)));

selfStencil=ones(P+1);
selfStencil=kron(eye(N),selfStencil);
upwindStencil=ones(P+1);
upwindStencil(1:2:P+1,:)=-1;
temporal = zeros(N);
temporal(2:N,1:N-1)=eye(N-1);
upwindStencil=kron(temporal, upwindStencil);
%Numerical flux matrix F, FU-initial condition is the flux
Flux=selfStencil + upwindStencil;

[currentUpper, nextLower, currentLower, previousUpper, coefficients] = centralFluxMatrices(N,P);
boundaryCondition=zeros(N*(P+1),1);
boundaryCondition(1:P+1)=(-1).^(0:P);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Expclicit time integration using 4th order Runge-Kutta%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t=0:0.001:6;
startT=t(1);
endT=t(end);
nT=length(t);
disp(nT);
deltaT=(endT-startT)/(nT-1);
results=zeros(N*(P+1),nT);
Results=zeros(N*(P+1),nT);
RESULTS=zeros(N*(P+1),nT);
IC=reshape(ApprLP(@initialCondition,X,P),N*(P+1),1);
Memory=zeros(N*(P+1),4);
MEM=zeros(N*(P+1),5);
results(:,1)=IC;
Results(:,1)=IC;
RESULTS(:,1)=IC;
tic
for i = 1:nT
    
    U=results(:,i);
    
    %Runge-kutta RK4
    %y(n+1)=y(n)+1/6*(k1+2*k2+2*k3+k4)
    %k1 = h*f( y(n) )
    %k2 = h*f( y(n)+k1/2 )
    %k3 = h*f( y(n)+k2/2 )
    %k4 = h*f( y(n)+k3 )
    %f(U)=source(U) +  -flux(U)
    
     
%      k1 = deltaT.*(iM*generalIntegral(N,U,@S,qW,qP,LP,LP) + ...
%          iM*iJ*generalIntegral(N,U,@F,qW,qP,LP,LPder) ...
%          - iM*iJ*(Flux*U - bc(i,deltaT).*boundaryCondition));
%      
%      k2 = deltaT.*(iM*generalIntegral(N,U+k1/2,@S,qW,qP,LP,LP) ...
%          +  iM*iJ*generalIntegral(N,U+k1/2,@F,qW,qP,LP,LPder) ...
%          - iM*iJ*(Flux*(U+k1/2) - bc(i,deltaT).*boundaryCondition ));
%      
%      k3 = deltaT.*(iM*generalIntegral(N,U+k2/2,@S,qW,qP,LP,LP) ...
%          + iM*iJ* generalIntegral(N,U+k2/2,@F,qW,qP,LP,LPder) ...
%          - iM*iJ*(Flux*(U+k2/2) - bc(i,deltaT).*boundaryCondition ));
%      
%      k4 = deltaT.*(iM*generalIntegral(N,U+k3,@S,qW,qP,LP,LP) ...
%          +  iM*iJ*generalIntegral(N,U+k3,@F,qW,qP,LP,LPder) ...
%          - iM*iJ*(Flux*(U+k3) - bc(i,deltaT).*boundaryCondition ));



%      k1 = deltaT.*(iM*generalIntegral(N,U,@S,qW,qP,LP,LP) + ...
%          iM*iJ*generalIntegral(N,U,@F,qW,qP,LP,LPder) ...
%          - iM*iJ*(upwindFlux(N,U,@F) - F(bc(i,deltaT)).*boundaryCondition));
%      
%      k2 = deltaT.*(iM*generalIntegral(N,U+k1/2,@S,qW,qP,LP,LP) ...
%          +  iM*iJ*generalIntegral(N,U+k1/2,@F,qW,qP,LP,LPder) ...
%          - iM*iJ*(upwindFlux(N,U+k1/2,@F) - F(bc(i,deltaT)).*boundaryCondition));
%      
%      k3 = deltaT.*(iM*generalIntegral(N,U+k2/2,@S,qW,qP,LP,LP) ...
%          + iM*iJ* generalIntegral(N,U+k2/2,@F,qW,qP,LP,LPder) ...
%          - iM*iJ*(upwindFlux(N,U+k2/2,@F) - F(bc(i,deltaT)).*boundaryCondition));
%      
%      k4 = deltaT.*(iM*generalIntegral(N,U+k3,@S,qW,qP,LP,LP) ...
%          +  iM*iJ*generalIntegral(N,U+k3,@F,qW,qP,LP,LPder) ...
%          - iM*iJ*(upwindFlux(N,U+k3,@F) - F(bc(i,deltaT)).*boundaryCondition));


%      k1 = deltaT.*(iM*generalIntegral(N,U,@S,qW,qP,LP,LP) + ...
%          iM*iJ*generalIntegral(N,U,@F,qW,qP,LP,LPder) ...
%          - iM*iJ*(LaxFriedrichsFlux(N,U,@F,@dF,bc(i,deltaT), currentUpper, nextLower, currentLower, previousUpper)));
%      
%      k2 = deltaT.*(iM*generalIntegral(N,U+k1/2,@S,qW,qP,LP,LP) ...
%          +  iM*iJ*generalIntegral(N,U+k1/2,@F,qW,qP,LP,LPder) ...
%          - iM*iJ*(LaxFriedrichsFlux(N,U+k1/2,@F,@dF,bc(i,deltaT), currentUpper, nextLower, currentLower, previousUpper)));
%      
%      k3 = deltaT.*(iM*generalIntegral(N,U+k2/2,@S,qW,qP,LP,LP) ...
%          + iM*iJ* generalIntegral(N,U+k2/2,@F,qW,qP,LP,LPder) ...
%          - iM*iJ*(LaxFriedrichsFlux(N,U+k2/2,@F,@dF,bc(i,deltaT), currentUpper, nextLower, currentLower, previousUpper)));
%      
%      k4 = deltaT.*(iM*generalIntegral(N,U+k3,@S,qW,qP,LP,LP) ...
%          +  iM*iJ*generalIntegral(N,U+k3,@F,qW,qP,LP,LPder) ...
%          - iM*iJ*(LaxFriedrichsFlux(N,U+k3,@F,@dF,bc(i,deltaT), currentUpper, nextLower, currentLower, previousUpper)));


     k1 = deltaT.*(iM*generalIntegral(N,U,@S,qW,qP,LP,LP) + ...
         iM*iJ*generalIntegral(N,U,@F,qW,qP,LP,LPder) ...
         - iM*iJ*(GodunovFlux(N,U,@F,bc(i,deltaT), currentUpper, nextLower, currentLower, previousUpper)));
     
     k2 = deltaT.*(iM*generalIntegral(N,U+k1/2,@S,qW,qP,LP,LP) ...
         +  iM*iJ*generalIntegral(N,U+k1/2,@F,qW,qP,LP,LPder) ...
         - iM*iJ*(GodunovFlux(N,U+k1/2,@F,bc(i,deltaT), currentUpper, nextLower, currentLower, previousUpper)));
     
     k3 = deltaT.*(iM*generalIntegral(N,U+k2/2,@S,qW,qP,LP,LP) ...
         + iM*iJ* generalIntegral(N,U+k2/2,@F,qW,qP,LP,LPder) ...
         - iM*iJ*(GodunovFlux(N,U+k2/2,@F,bc(i,deltaT), currentUpper, nextLower, currentLower, previousUpper)));
     
     k4 = deltaT.*(iM*generalIntegral(N,U+k3,@S,qW,qP,LP,LP) ...
         +  iM*iJ*generalIntegral(N,U+k3,@F,qW,qP,LP,LPder) ...
         - iM*iJ*(GodunovFlux(N,U+k3,@F,bc(i,deltaT), currentUpper, nextLower, currentLower, previousUpper)));
     
     results(:,i+1) = U +(1/6)*(k1 + 2*k2 + 2*k3 + k4);
    
     
    if i/floor(nT/10)==floor(i/floor(nT/10))
        fprintf('%s','*');
    end    
end
fprintf('%s\n', '');
disp("Runge-Kutta: " + toc)

tic
for i=1:nT
     W=Results(:,i);
%     Results(:,i+1)= W + deltaT*Mass\Jacobian\(Jacobian*generalIntegral(N,W,@S,qW,qP,LP,LP) ...
%         +  generalIntegral(N,W,@F,qW,qP,LP,LPder) ...
%         - Flux*W + bc(i,deltaT).*boundaryCondition);
    
%     Results(:,i+1)=W + deltaT.*(iM*generalIntegral(N,W,@S,qW,qP,LP,LP) + ...
%           iM*iJ*generalIntegral(N,W,@F,qW,qP,LP,LPder) ...
%           - iM*iJ*(Flux*W - bc(i,deltaT).*boundaryCondition));
      
      %LaxFriedrichsFlux(N,U,F,dF,bc,currentUpper, nextLower, currentLower, previousUpper)
%       Results(:,i+1)=W + deltaT.*(iM*generalIntegral(N,W,@S,qW,qP,LP,LP) + ...
%           iM*iJ*generalIntegral(N,W,@F,qW,qP,LP,LPder) ...
%           - iM*iJ*(LaxFriedrichsFlux(N,W,@F,@dF,bc(i,deltaT),currentUpper, nextLower, currentLower, previousUpper)));

%     Results(:,i+1)=W + deltaT.*(iM*generalIntegral(N,W,@S,qW,qP,LP,LP) + ...
%           iM*iJ*generalIntegral(N,W,@F,qW,qP,LP,LPder) ...
%           - iM*iJ*(upwindFlux(N,W,@F) - F(bc(i,deltaT)).*boundaryCondition));

%     Results(:,i+1)=W + deltaT.*(iM*generalIntegral(N,W,@S,qW,qP,LP,LP) + ...
%           iM*iJ*generalIntegral(N,W,@F,qW,qP,LP,LPder) ...
%           - iM*iJ*(flux(N,W,@F) - bc(i,deltaT).*boundaryCondition));


      Results(:,i+1)=W + deltaT.*(iM*generalIntegral(N,W,@S,qW,qP,LP,LP) + ...
          iM*iJ*generalIntegral(N,W,@F,qW,qP,LP,LPder) ...
          - iM*iJ*(GodunovFlux(N,W,@F,bc(i,deltaT), currentUpper, nextLower, currentLower, previousUpper)...
          ));
    
    if i/floor(nT/10)==floor(i/floor(nT/10))
        fprintf('%s','*');
    end 
end    
fprintf('%s\n', '');
disp("Euler: " + toc)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%4th order Adams-Bashfort%%%%%%%%
tic
for i=1:3
    
    R=RESULTS(:,i);
     k1 = deltaT.*(iM*generalIntegral(N,R,@S,qW,qP,LP,LP) + ...
         iM*iJ*generalIntegral(N,R,@F,qW,qP,LP,LPder) ...
         - iM*iJ*(GodunovFlux(N,R,@F,bc(i,deltaT), currentUpper, nextLower, currentLower, previousUpper)));
     
     k2 = deltaT.*(iM*generalIntegral(N,R+k1/2,@S,qW,qP,LP,LP) ...
         +  iM*iJ*generalIntegral(N,R+k1/2,@F,qW,qP,LP,LPder) ...
         - iM*iJ*(GodunovFlux(N,R+k1/2,@F,bc(i,deltaT), currentUpper, nextLower, currentLower, previousUpper)));
     
     k3 = deltaT.*(iM*generalIntegral(N,R+k2/2,@S,qW,qP,LP,LP) ...
         + iM*iJ* generalIntegral(N,R+k2/2,@F,qW,qP,LP,LPder) ...
         - iM*iJ*(GodunovFlux(N,R+k2/2,@F,bc(i,deltaT), currentUpper, nextLower, currentLower, previousUpper)));
     
     k4 = deltaT.*(iM*generalIntegral(N,R+k3,@S,qW,qP,LP,LP) ...
         +  iM*iJ*generalIntegral(N,R+k3,@F,qW,qP,LP,LPder) ...
         - iM*iJ*(GodunovFlux(N,R+k3,@F,bc(i,deltaT), currentUpper, nextLower, currentLower, previousUpper)));
     RESULTS(:,i+1) = R +(1/6)*(k1 + 2*k2 + 2*k3 + k4);
end    
% Memory(:,2)=RESULTS(:,1);
% Memory(:,3)=RESULTS(:,2);
% Memory(:,4)=RESULTS(:,3);

Memory(:,2)=deltaT.*(iM*generalIntegral(N,RESULTS(:,1),@S,qW,qP,LP,LP) + ...
    iM*iJ*generalIntegral(N,RESULTS(:,1),@F,qW,qP,LP,LPder) ...
    - iM*iJ*(GodunovFlux(N,RESULTS(:,1),@F,bc(i,deltaT), currentUpper, nextLower, currentLower, previousUpper)));
Memory(:,3)=deltaT.*(iM*generalIntegral(N,RESULTS(:,2),@S,qW,qP,LP,LP) + ...
    iM*iJ*generalIntegral(N,RESULTS(:,2),@F,qW,qP,LP,LPder) ...
    - iM*iJ*(GodunovFlux(N,RESULTS(:,2),@F,bc(i,deltaT), currentUpper, nextLower, currentLower, previousUpper)));
Memory(:,4)=deltaT.*(iM*generalIntegral(N,RESULTS(:,3),@S,qW,qP,LP,LP) + ...
    iM*iJ*generalIntegral(N,RESULTS(:,3),@F,qW,qP,LP,LPder) ...
    - iM*iJ*(GodunovFlux(N,RESULTS(:,3),@F,bc(i,deltaT), currentUpper, nextLower, currentLower, previousUpper)));

for i=4:nT
    R=RESULTS(:,i);
    Memory(:,1)=Memory(:,2);
    Memory(:,2)=Memory(:,3);
    Memory(:,3)=Memory(:,4);
    Memory(:,4)=deltaT.*(iM*generalIntegral(N,R,@S,qW,qP,LP,LP) + ...
         iM*iJ*generalIntegral(N,R,@F,qW,qP,LP,LPder) ...
         - iM*iJ*(GodunovFlux(N,R,@F,bc(i,deltaT), currentUpper, nextLower, currentLower, previousUpper)));
    RESULTS(:,i+1)=R+55/24*Memory(:,4)-59/24*Memory(:,3)+37/24*Memory(:,2)-9/24*Memory(:,1);
    if i/floor(nT/10)==floor(i/floor(nT/10))
        fprintf('%s','*');
    end 
end    
fprintf('%s\n', '');
disp("Adams-Bashfort: " + toc);

fps=30;
frameT=floor(nT/(fps*(endT-startT)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%Plotting%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:frameT:nT
    tic
    U=reshape(results(:,i),P+1,N);
    W=reshape(Results(:,i),P+1,N);
    RR=reshape(RESULTS(:,i),P+1,N);
    Y=plotLP(x,X,U);
    Z=plotLP(x,X,W);
    ABY=plotLP(x,X,RR);
    plot(x,Y,'r');
    hold on
    plot(x,Z,'b');
    plot(x,ABY,'g');
    hold off
    axis([a b -1.5 1.5]);
    text(1.02,b+0.1,'Time');
    text(1.02,1,num2str((i-1)*(endT-startT)*1/nT));
    %disp(1/fps-toc)
    %pauseT=max( 1/fps-toc,0 )
    pause(max( 1/fps-toc,0 ));
end    

% for i=1:10:nT
%     tic
%     U=reshape(results(:,i),P+1,N);
%     %W=reshape(Results(:,i),P+1,N);
%     Y=plotLP(x,X,U);
%     %Z=plotLP(x,X,W);
%     plot(x,Y,'r');
%     %hold on
%     %plot(x,Z,'b');
%     %hold off
%     axis([a b -1.5 1.5]);
%     text(1.02,b+0.1,'Time');
%     text(1.02,1,num2str((i-1)*(endT-startT)*1/nT));
%     %disp(1/fps-toc)
%     %pauseT=max( 1/fps-toc,0 )
%     pause(0.01);
% end    


% for i=2001:-200:1
%     U=reshape(results(:,i),P+1,N);
%     Y=plotLP(x,X,U);
%     W=reshape(Results(:,i),P+1,N);
%     Z=plotLP(x,X,W);
%     
%     figure
%     plot(x,Y,'r');
%     hold on
%     plot(x,Z,'b')
%     axis([a b -1.5 1.5]);
% end

function y=F(x)
    %y=x+0.1*x.^2;
    %y=5*x+1.5*x.^2;
    y=x;
    %y=x+0.5*exp(-3*x-5)-0.05*x.^2;
    %y=0.01*x.^2;
    %y=exp(0.2*x)+x;
    %y=nthroot(x,3);
    %y=exp(-x.^2);
    %y=0.5*abs(x);
   
end

function y=dF(x)
    y=1+0.2*x;
    %y=ones(size(x));
    %y=0.2*exp(0.2*x)+1;
    %y=5+2.8*x;
end

function y=S(x)
    y=-0.2*x;
    %y=-0.5*x.*abs(x);
    %y=x-x;
    %y=zeros(size(x));
    %y=x;
end

function y=bc(i,deltaT)
    %y=sin(10*deltaT*(i-1));
    y=sin(4*deltaT*(i-1));
    %y=floor(i/100)-2*floor(i/200);
    %y=0.1;
end  

function y=initialCondition(x)
    %y=sin(3*x);
    y=-0.5*exp(-10*(x-0.2).^2);
end


















