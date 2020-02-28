function y=fluxTerm(N,A,U,F,bcA,bcU,ecA,ecU,Beta,a0,rho, currentUpper, nextLower, currentLower, previousUpper)
    P=length(U)/N-1;

    Ar=(-1).^(0:P)*A(1:P+1);%Value of A at the inlet
    Ur=(-1).^(0:P)*U(1:P+1);%value of u at the inlet
    
%     upperBound=flux(F,currentUpper*A,currentUpper*U,nextLower*A+[zeros(N-1,1);ecA],...
%         nextLower*U+0*[zeros(N-1,1);ecU],Beta,a0,rho)-F(currentUpper*A,currentUpper*U);
%     lowerBound=flux(F,previousUpper*A+Aboundary(bcA,Ar)*[1;zeros(N-1,1)],previousUpper*U+Uboundary(bcU,Ur)*[1;zeros(N-1,1)], ...
%         currentLower*A,currentLower*U,Beta,a0,rho)-F(currentLower*A,currentLower*U);

    upperBound=flux(F,currentUpper*A,currentUpper*U,nextLower*A+[zeros(N-1,1);ecA],...
        nextLower*U+1*[zeros(N-1,1);ecU],Beta,a0,rho)-F(currentUpper*A,currentUpper*U);
    lowerBound=flux(F,previousUpper*A+bcA*[1;zeros(N-1,1)],previousUpper*U+bcU*[1;zeros(N-1,1)], ...
        currentLower*A,currentLower*U,Beta,a0,rho)-F(currentLower*A,currentLower*U);

    y=kron(upperBound,ones(P+1,1))...
        -kron(ones(N,1), (-1).^(0:P)').*...
        kron(lowerBound,ones(P+1,1));
    
end



function y=flux(F,a1,u1,a2,u2,Beta,a0,rho)
    y=F(Aupwinded(a1,u1,a2,u2,Beta,a0,rho),Uupwinded(a1,u1,a2,u2,Beta,a0,rho));
end



function y=W1(a,u,Beta,a0,rho)
    %u0=1;
    u0=0;
    gamma=sqrt(Beta./(2*rho*a0));
    %y=u-u0+4*(gamma.*(nthroot(a,4) - nthroot(a0,4)));
    
    %Faster, but can output imaginary roots:
    y=u-u0+4*(gamma.*(nthroot(a,4) - a0.^0.25));
    
    %Computational modelling of 1D blood...
%     gamma=sqrt(Beta./(2*rho));
%     y=u-u0+4*gamma.*(nthroot(a,4));
end

function y=W2(a,u,Beta,a0,rho)
    %u0=1;
    u0=0;
    gamma=sqrt(Beta./(2*rho*a0));
    %y=u-u0-4*(gamma.*(nthroot(a,4) - nthroot(a0,4)));
    
    %Faster, but can output imaginary roots:
    y=u-u0-4*(gamma.*(a.^0.25 - a0.^0.25));
    
    %Computational modelling of 1D blood...
%     gamma=sqrt(Beta./(2*rho));
%     y=u-u0-4*gamma.*(nthroot(a,4));
end

function y=Aupwinded(a1,u1,a2,u2,Beta,a0,rho)
    %y=nthroot( (sqrt(rho.*a0./(Beta*32)).*(W1(a1,u1,Beta,a0,rho)-W2(a2,u2,Beta,a0,rho)) + nthroot(a0,4))  ,4);
    %y=(sqrt(rho.*a0./(Beta*32)).*(W1(a1,u1,Beta,a0,rho)-W2(a2,u2,Beta,a0,rho))+nthroot(a0,4)).^4;
    
    
    %Faster, but can output imaginary roots:
    y=(sqrt(rho.*a0./(Beta*32)).*(W1(a1,u1,Beta,a0,rho)-W2(a2,u2,Beta,a0,rho))+a0.^0.25).^4;
    
    %Computational modelling of 1D blood...
    %y=(sqrt(rho./(Beta)).*(W1(a1,u1,Beta,a0,rho)-W2(a2,u2,Beta,a0,rho))).^4;
end

function y=Uupwinded(a1,u1,a2,u2,Beta,a0,rho)
%Added +U0? to the end
    y=0.5*(W1(a1,u1,Beta,a0,rho)+W2(a2,u2,Beta,a0,rho));%+1;

end

function Al=Aboundary(Abc,Ar)
    %Al=(2*nthroot(Abc,4)-nthroot(Ar,4)).^4;
    Al=(2*Abc.^0.25-Ar.^0.25).^4;
end

function Ul=Uboundary(Ubc,Ur)
    Ul=2*Ubc-Ur;
end