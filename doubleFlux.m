function y=doubleFlux(N,U,F,bc, currentUpper, nextLower, currentLower, previousUpper)
    P=length(U)/N-1;
    upperBound=GoFlux(F,currentUpper*U,nextLower*U)-F(currentUpper*U);
    lowerBound=GoFlux(F,previousUpper*U+bc*[1 zeros(1,N-1)]',currentLower*U)-F(currentLower*U);
    y=kron(upperBound,ones(P+1,1))...
        -kron(ones(N,1), (-1).^(0:P)').*...
        kron(lowerBound,ones(P+1,1));
end