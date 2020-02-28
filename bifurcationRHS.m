function [RHSA, RHSU]=bifurcationRHS(N,A,U,qW,qP,LPminus,LPplus,LP,epsilon,boundaries,bifurcations,parameters, currentUpper, nextLower, currentLower, previousUpper,iM,iJ)
    %
    %INPUTS:
    %
    %N - the number of elements
    %
    %boundaryConditions
    %parameters(i,:)=[Beta a0 rho]
    
    segments=size(parameters,1);
    P=length(U)/N-1;
    B=size(bifurcations,1);
    rho=parameters(1,3);
    pext=parameters(1,7);

    for i=1:segments
        Atemp=A(i,:)';
        Utemp=U(i,:)';
%         disp("A: " + size(Atemp))
%         disp("U: " + size(Utemp))
        %Notice that returns all segments as row vectors instead of
        %columns. Should be compatible with rank 3 resultsA and resultsU
        RHSA(i,:)=RHS2(N,Atemp,Utemp,@F1,@S1,qW,qP,LPminus,LPplus,LP,epsilon,boundaries(i,:),parameters(i,:),...
            currentUpper, nextLower, currentLower, previousUpper,iM,iJ)';
        RHSU(i,:)=RHS2(N,Atemp,Utemp,@(a,u) F2(a,u,pext,parameters(i,1),parameters(i,2),rho),...
        @(a,u) S2(a,u,rho),qW,qP,LPminus,LPplus,LP,epsilon,boundaries(i,:),parameters(i,:),...
            currentUpper, nextLower, currentLower, previousUpper,iM,iJ)';
    end    
end


function y=RHS2(N,A,U,F,S,qW,qP,LPminus,LPplus,LP,epsilon,boundaries,parameters, currentUpper, nextLower, currentLower, previousUpper,iM,iJ)
    %
    %INPUTS:
    %
    %N - the number of elements
    y=-iM*iJ*derivativeIntegral2(N,A,U,F,qW,qP,LPminus,LPplus,LP,epsilon) ...
        -iJ*iM*fluxTerm(N,A,U,F,boundaries(1),boundaries(2),boundaries(3),boundaries(4),parameters(1),parameters(2),parameters(3), currentUpper, nextLower, currentLower, previousUpper)...
        +iM*generalIntegral2(N,A,U,S,qW,qP,LP,LP);
end

function y=F1(a,u)
    y=a.*u;
end

function y=F2(a,u,pext,Beta,a0,rho)
    y=u.*u/2+(pext+(Beta./a0).*(sqrt(a)-sqrt(a0)))/rho;
    
    %Computational modelling of 1D blood...
    %y=u.*u/2+(pext+(Beta)./a0.*(sqrt(a)-sqrt(a0)))/rho;
end  

function y=S1(a,u)
    y=zeros(size(a));
end

function y=S2(a,u,rho)
    y=-22*pi*0.000025*u ./(rho*a);
    %y=-22*pi*u ./(rho*a);
    %y=zeros(size(a));
end