close all;
clear all;
a=0;
b=1;
%resolution of the plot, solution is evaluated at resolution points
resolution = 1000;
h=(b-a)/resolution;
x=a:h:b;

N=4;
P=6;
deltaX=(b-a)/N;
X=a:deltaX:b;

%Solve advection-equation du/dt+cdu/dx=s(u), x \in [a,b]
%u(0,t)=sin(t), u(x,0)=0

c=1;


    %4th order Gauss-Legendre
    qW=[0.3478548451374538; 0.6521451548625461; 0.6521451548625461; 0.3478548451374538];
    qP=[-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526];
    LP=zeros(11,4);
    LP(1,:)=1;
    LP(2,:)=qP;
    LP(3,:)=1/2*(3*qP.^2-1);
    LP(4,:)=1/2*(5*qP.^3-3.*qP);
    LP(5,:)=1/8*(35*qP.^4-30*qP.^2+3);
    LP(6,:)=1/8*(63*qP.^5-70*qP.^3+15*qP);
    LP(7,:)=1/16*(231*qP.^6-315*qP.^4 ...
        +105*qP.^2-5);
    
    LP(8,:)=1/16*(429*qP.^7-693*qP.^5 ...
        +315*qP.^3-35*qP);
   
    LP(9,:)=1/128*(6435*qP.^8-12012*qP.^6 ...
        +6930*qP.^4-1260*qP.^2+35);
    
    LP(10,:)=1/128*(12155*qP.^9-25740*qP.^7 ...
        +18018*qP.^5-4620*qP.^3+315*qP);
   
    LP(11,:)=1/256*(46189*qP.^10-109395*qP.^8 ...
        +90090*qP.^6-30030*qP.^4+3465*qP.^2-63);

%chaotic settings: a,b=0,20 N=32,P=5 nT=1000


selfStencil=toeplitz((-1).^(0:P),ones(1,P+1));
selfStencil=kron(eye(N),selfStencil);
upwindStencil=ones(P+1);
upwindStencil(1:2:P+1,:)=-1;
temporal = zeros(N);
temporal(2:N,1:N-1)=eye(N-1);
upwindStencil=kron(temporal, upwindStencil);
%Numerical flux matrix F, FU-initial condition is the flux
F=selfStencil + upwindStencil;%F-S

massMatrix=diag(2./(2.*(0:P)+1));
massMatrix=kron(diag(diff(X)/2),massMatrix);

U=zeros((P+1)*N,1);

boundaryCondition=zeros(N*(P+1),1);
boundaryCondition(1:P+1)=(-1).^(0:P);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Explicit time stepping using Ringe-Kutta%%%%%%%%%%%%%%%%%%%%%%%%%
startT=0;
endT=5;
nT=floor((endT-startT)/0.001)
%nT=5000;%10000
deltaT=(endT-startT)/(nT-1);
t=startT:deltaT:endT;

results=zeros(N*(P+1),nT);
Results=zeros(N*(P+1),nT);
%results(1:P+1:end,1)=1;
%results(:,1)=reshape(ApprLP(@sin,x,X,P),N*(P+1),1);
tic
for i = 1:nT
    %results(:,i)=results(:,i-1) + deltaT*massMatrix\(-c*(F)*results(:,i-1)+c*sin((i-1)*deltaT).*boundaryCondition)
    %results(:,i)=results(:,i-1) + deltaT*massMatrix\(-F*results(:,i-1));
    %results(:,i+1)=results(:,i) + deltaT*massMatrix\(-F*results(:,i));
    U=results(:,i);
    
    
    %Source term
    Q=zeros(N*(P+1),1);
%     for e=1:N;
%         Q(1+(e-1)*(P+1):e*(P+1))=deltaX/2.*sourceIntegral(U(1+(e-1)*(P+1):e*(P+1)),@S);
%     end
    %Q=sourceIntegral2(diff(X),U,@S, qW, qP, LP);
    Q=diag(kron(diff(X)./2,ones(1,P+1)))*(generalIntegral(N,U,@S,qW,qP,LP,LP));
    
    k1 = massMatrix\(Q-c*F*U + c*bc(i,deltaT).*boundaryCondition ).*deltaT;
    k2 = massMatrix\(Q-c*F* (U + k1 /2) + c*bc(i,deltaT).*boundaryCondition).*deltaT;
    k3 = massMatrix\(Q-c*F* (U + k2 /2) + c*bc(i,deltaT).*boundaryCondition).*deltaT;
    k4 = massMatrix\(Q-c*F* (U + k3) + c*bc(i,deltaT).*boundaryCondition).*deltaT;
    results(:,i+1) = U +(1/6)*(k1 + 2*k2 + 2*k3 + k4);
    
end
time=toc;
fprintf('%s  %f\n', "Runge-Kutta: " + time);

tic

for i = 1:nT
    W=Results(:,i);
    Results(:,i+1)= W + deltaT*(massMatrix\(Q -c*F*W + c*bc(i,deltaT).*boundaryCondition ));
    
end
time = toc;
disp(" ")
fprintf('%s  %f\n', "Euler: " + time);
disp(" ")



fps=30;
frameT=floor(nT/(fps*(endT-startT)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:frameT:nT
    tic
    U=reshape(results(:,i),P+1,N);
    W=reshape(Results(:,i),P+1,N);
    Y=plotLP(x,X,U);
    %Z=plotLP(x,X,W);
    plot(x,Y,'r');
    hold on
    %plot(x,Z,'b');
    hold off
    axis([a b -1.5 1.5]);
    text(1.02,b+0.1,'Time')
    text(1.02,1,num2str((i-1)*(endT-startT)*1/(nT-1)));
    %disp(1/fps-toc)
    %pauseT=max( 1/fps-toc,0 )
    pause(max( 1/fps-toc,0 ));
end    

function y=bc(i,deltaT)
    %y=cos(10*deltaT*(i-1));
    y=sin(10*deltaT*(i-1));
    %y=floor(i/100)-2*floor(i/200);
end    

function y=S(x)
    %y=-0.2*x.*abs(x);
    %y=-x+1;
    y=-x;
    
end    
    













