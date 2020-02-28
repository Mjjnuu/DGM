function int = sourceIntegral2(deltaX,U,S, qW,qP,LP)
    %U is vector of coefficients in Legendre basis in one element
    %S is the source function
    %Calculates integrals \int_I S(u_h)L_q dx, where I=[-1,1], 
    %u_h=\sum_{p=0}^P U(p)L_p, and L_i is the Legendre polynomial
    %of degree i \in {0,1,...,P}.
    N=length(deltaX);
    P=length(U)/N-1;
      
    %size(diag(kron(deltaX./2,ones(1,P+1))))                    % 50 50
    %size(reshape(U,P+1,N)'*LP(1:P+1,:))                        % 10 4
    %size(S(reshape(U,P+1,N)'*LP(1:P+1,:))*diag(qW))            % 10 4
    %size((kron(S(reshape(U,P+1,N)'*LP(1:P+1,:))*diag(qW),ones(P+1,1)).*kron(ones(N,1),LP(1:P+1,:)))) %50 4
    %size(ones(N*(P+1),1)) %50 1
    
    int=diag(kron(deltaX./2,ones(1,P+1)))* ...
        ((kron(S(reshape(U,P+1,N)'*LP(1:P+1,:))*diag(qW),ones(P+1,1)).*kron(ones(N,1),LP(1:P+1,:)))*ones(length(qW),1));
       
end