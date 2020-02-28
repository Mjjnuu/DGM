function y=centralFlux(N,U,F)

    P=length(U)/N-1;
    [A,B,C,D,E]=centralFluxMatrices(N,P);
    %y=1/2*(F(A*U)+F(B*U)-E*(F(C*U)+F(D*U)));
    y=1/2*((A+B-E.*(C+D))*U);

end