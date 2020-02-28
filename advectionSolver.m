function results = advectionSolver(t,x,X,c,P,bc)
    %Input:
    %t, the time domain vector. Must satisfy t(i)-t(i-1)=t(j)-t(j-1)
    %for all i and j.
    %x, the spatial domain vector
    %X, partition of spatial domain
    %c, speed
    %P, maximum degree of polynomial space
    %bc, vector containing boundary conditions u(x(1),t)
    %length(bc)=length(t)

    %Returns:
    %
    
    %Solves the advection equation du/dt+c*du/dx=0 using discontinuous
    %Galerkin method and 4th order Runge-Kutta time integration.
    
    a=x(1);
    b=x(end);
    
    N=length(X)-1;
    startT=t(1);
    endT=t(end);
    nT=length(t);
    deltaT=(endT-startT)/nT;
    
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
    
    
    %Explicit time integration using 4th order Runge-Kutta method:
    results=zeros(N*(P+1),nT);
    
    for i = 1:nT
        U=results(:,i);
        k1 = massMatrix\(-c*F*U + c*bc(i).*boundaryCondition ).*deltaT;
        k2 = massMatrix\(-c*F* (U + k1 /2) + c*bc(i).*boundaryCondition).*deltaT;
        k3 = massMatrix\(-c*F* (U + k2 /2) + c*bc(i).*boundaryCondition).*deltaT;
        k4 = massMatrix\(-c*F* (U + k3) + c*bc(i).*boundaryCondition).*deltaT;
        results(:,i+1) = U +(1/6)*(k1 + 2*k2 + 2*k3 + k4);
        
    end
    


end