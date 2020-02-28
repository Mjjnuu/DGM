function y=LaxFriedrichsFlux(N,U,F,dF,bc, currentUpper, nextLower, currentLower, previousUpper)
    P=length(U)/N-1;
    %disp(size());
    y=kron(LFflux(F,dF,nextLower*U,currentUpper*U),ones(P+1,1))...
        -kron(ones(N,1), (-1).^(0:P)').*...
        kron(LFflux(F,dF,currentLower*U,previousUpper*U+bc*[1 zeros(1,N-1)]'),ones(P+1,1));
        
end