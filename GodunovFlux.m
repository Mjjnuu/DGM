function y=GodunovFlux(N,U,F,bc, currentUpper, nextLower, currentLower, previousUpper)
    P=length(U)/N-1;
    y=kron(GoFlux(F,currentUpper*U,nextLower*U),ones(P+1,1))...
        -kron(ones(N,1), (-1).^(0:P)').*...
        kron(GoFlux(F,previousUpper*U+bc*[1 zeros(1,N-1)]',currentLower*U),ones(P+1,1));


end

