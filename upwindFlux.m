function y=upwindFlux(N,U,F)
    P=length(U)/N-1;
    aaa=kron(diag(ones(1,N)),ones(P+1));
    
    bbb=kron(toeplitz([0 1 zeros(1,N-2)],zeros(1,N)), ones(P+1));
    ddd=diag(kron(ones(1,N), (-1).^(0:P)));
    
    y=F(aaa*U)-ddd*F(bbb*U);

end