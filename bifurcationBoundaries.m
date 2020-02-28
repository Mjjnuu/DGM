function newBoundaries=bifurcationBoundaries(N,A,U,boundaries,bifurcations,parameters,boundaryConditions)
    newBoundaries=boundaries;
    segments=size(parameters,1);
    P=length(U)/N-1;
    B=size(bifurcations,1);
    rho=parameters(3,1);
    Q=zeros(1,6);
    
    %loop over all bifurcations
    for i=1:B
        %segments taking part in the bifurcation
        seg=bifurcations(i,:);
        %A and u at the end of the segment seg(1):
        Q(1)=ones(1,P+1)*A(seg(1),(N-1)*(P+1)+1:N*(P+1))';
        Q(2)=ones(1,P+1)*U(seg(1),(N-1)*(P+1)+1:N*(P+1))';
        %A and u at the beginning of the segment seg(2):
        Q(3)=(-1).^(0:P)*A(seg(2),1:P+1)';
        Q(4)=(-1).^(0:P)*U(seg(2),1:P+1)';
        %A and u at the beginning of the segment seg(3):
        Q(5)=((-1).^(0:P))*A(seg(3),1:P+1)';
        Q(6)=((-1).^(0:P))*U(seg(3),1:P+1)';
        
        %Upwinded A1,u1,A2,u2,A3,u3:
        temp=FastNewtonRhapson(Q,10,Q,parameters(seg,2),parameters(seg,1),rho);
        newBoundaries(seg(1),3:4)=temp(1:2);
        newBoundaries(seg(2),1:2)=temp(3:4);
        newBoundaries(seg(3),1:2)=temp(5:6);
    end    
    
    %
    Ar=(-1).^(0:P)*A(1,1:P+1)';
    Ur=(-1).^(0:P)*U(1,1:P+1)';
    newBoundaries(1,1)=Aboundary(boundaryConditions(1),Ar);
    newBoundaries(1,2)=Uboundary(boundaryConditions(2),Ur);
%     newBoundaries(1,1)=boundaryConditions(1);
%     newBoundaries(1,2)=boundaryConditions(2);
    
    %Implementation of terminal conditions here
    
end

%From fluxTerm.m
%Calculates what the value of A on the upper end of the "virtual" segment
%left from the first segment should be so that the upwinded value Aupwinded
%would be equal to the desired inflow boundary condition (coming from a file/measurement/user)
function Al=Aboundary(Abc,Ar)
    Al=(2*Abc.^0.25-Ar.^0.25).^4;
end
%Calculates what the value of U on the upper end of the "virtual" segment
%left from the first segment should be so that the upwinded value Uupwinded
%would be equal to the desired inflow boundary condition (coming from a file/measurement/user)
function Ul=Uboundary(Ubc,Ur)
    Ul=2*Ubc-Ur;
end

%Newton-Rhapson iteration for solving a system of equations determining the
%upwinded parameters. See https://spiral.imperial.ac.uk/handle/10044/1/8784
%equations 2.86-2.91
function y=FastNewtonRhapson(x0,iterations,Q,A0,Beta,rho)
    for i=1:iterations
        x1=x0'-FastNumJacob(x0,A0,Beta,rho)\fastF(x0,Q,A0,Beta,rho);
        x0=x1';
    end 
%     x1=x0'-FastNumJacob(x0,Q,A0,Beta,rho)\fastF(x0,Q,A0,Beta,rho);
%     x0=x1';
%     x1=x0'-FastNumJacob(x0,Q,A0,Beta,rho)\fastF(x0,Q,A0,Beta,rho);
%     x0=x1';
    
    y=x0;
end
%Exact Jacobian of the system
function Jf=FastNumJacob(x,A0,Beta,rho)
    Jf=zeros(6);
    Jf(1,1)=-sqrt(Beta(1)/(2*rho*A0(1)))*x(1).^(-0.75);
    Jf(1,2)=-1;
    
    Jf(2,3)=sqrt(Beta(2)/(2*rho*A0(2)))*x(3).^(-0.75);
    Jf(2,4)=-1;
    
    Jf(3,5)=sqrt(Beta(3)/(2*rho*A0(3)))*x(5).^(-0.75);
    Jf(3,6)=-1;
    
    Jf(4,1)=x(2);
    Jf(4,2)=x(1);
    Jf(4,3)=-x(4);
    Jf(4,4)=-x(3);
    Jf(4,5)=-x(6);
    Jf(4,6)=-x(5);
    
    Jf(5,1)=0.5*Beta(1)/A0(1)*x(1).^(-0.5);
    Jf(5,2)=rho*x(2);
    Jf(5,3)=-0.5*Beta(2)/A0(2)*x(3).^(-0.5);
    Jf(5,4)=-rho*x(4);
    
    Jf(6,1)=0.5*Beta(1)/A0(1)*x(1).^(-0.5);
    Jf(6,2)=rho*x(2);
    Jf(6,5)=-0.5*Beta(3)/A0(3)*x(5).^(-0.5);
    Jf(6,6)=-rho*x(6);
end

function y=fastF(x,Q,A0,Beta,rho)
    y(1,1)=Q(2)-x(2)+4*sqrt(Beta(1)/(2*rho*A0(1)))*(Q(1).^0.25-x(1).^0.25);
    y(2,1)=Q(4)-x(4)-4*sqrt(Beta(2)/(2*rho*A0(2)))*(Q(3).^0.25-x(3).^0.25);
    y(3,1)=Q(6)-x(6)-4*sqrt(Beta(3)/(2*rho*A0(3)))*(Q(5).^0.25-x(5).^0.25);
    y(4,1)=x(1).*x(2)-x(3).*x(4)-x(5).*x(6);
    y(5,1)=rho*0.5*(x(2).^2-x(4).^2)+Beta(1)/A0(1)*(sqrt(x(1))-sqrt(A0(1)))-Beta(2)/A0(2)*(sqrt(x(3))-sqrt(A0(2)));
    y(6,1)=rho*0.5*(x(2).^2-x(6).^2)+Beta(1)/A0(1)*(sqrt(x(1))-sqrt(A0(1)))-Beta(3)/A0(3)*(sqrt(x(5))-sqrt(A0(3)));
    
end
