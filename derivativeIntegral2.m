function int=derivativeIntegral2(N,A,U,F,qW,qP,LPminus,LPplus,V,epsilon)
    %Same as derivativeIntegral, but assumes that F takes two inputs
    P=length(U)/N-1;
    symmetricDifference=F(reshape(A,P+1,N)'*LPplus(1:P+1,:), reshape(U,P+1,N)'*LPplus(1:P+1,:))...
        -F(reshape(A,P+1,N)'*LPminus(1:P+1,:), reshape(U,P+1,N)'*LPminus(1:P+1,:));
    symmetricDifference=symmetricDifference/(2*epsilon);
    
    int=(kron(symmetricDifference*diag(qW),ones(P+1,1)).*kron(ones(N,1),V(1:P+1,:)))*ones(length(qP),1);

end