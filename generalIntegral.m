function int = generalIntegral(N,U,F,qW,qP,LP,V)
    %INPUTS:
    
    %deltaX is the vector containing the lengths of the elements.
    %size(deltaX)= [1,N].
    
    %U is vector of coefficients in Legendre basis 
    %U=[U_1;U_2;...;U_N], where U_i is the coefficient vector on element i.
    
    %F(U(x)) is a function
    
    %qW is a row vector containing quadrature weigths
    %qP is a row vector containing quadrature points
    %|qW|=|qP|=q
    
    %LP is a (P+1) x q matrix. LP(p,i) is the value of the Legendre 
    %polynomial of degree p-1 evaluated at i:th quadrature point i.e. 
    %LP(p,i)=L_(p-1)(qP(i)), where L_i is the
    %Legendre polynomial of degree i-1.
    
    %V is a (P+1) x q matrix. v_i , i \in {0,1,...,P} are the testing
    %functions, degree(v_i)=i. V(p+1,i) is the value of p+1:th testing
    %function evaluated at i:th quadrature point, i.e. V(p+1,i)=v_p(qP(i)).
    
  
    %PURPOSE:
    %Calculates the integrals \int_I F(u_h(x))*v_j(x) dx using numerical 
    %quadrature determined by qW and qP. Here u_h=\sum_{p=0}^P U_ep*L_p on 
    %element e, U_e is the coefficient vector with respect to Legendre 
    %polynomial basis on element e.  
    
    %int is a N*(P+1)x1 matrix. int=[I_1;I_2;...;I_N], where I_e are the
    %integral vectors of element e. I_e(j)=\int_e F(u_h(x))*v_j(x) dx =
    %\sum_{i=1}^q qW(i)*F(u_h(qP(i)))*v_j(qP(i))
    
    %N=length(deltaX);
    P=length(U)/N-1;
    
%     int=diag(kron(deltaX./2,ones(1,P+1)))* ...
%         ((kron(F(reshape(U,P+1,N)'*LP(1:P+1,:))*diag(qW),ones(P+1,1)).*kron(ones(N,1),LP(1:P+1,:)))*ones(length(qW),1));
    

    %int=(kron(F(reshape(U,P+1,N)'*LP(1:P+1,:))*diag(qW),ones(P+1,1)).*kron(ones(N,1),V(1:P+1,:)))*ones(length(qW),1);
    
    AA=F(reshape(U,P+1,N)'*LP(1:P+1,:));
    Kron=kron(AA*diag(qW),ones(P+1,1));
    size(AA);
    size(Kron);
    size(kron(ones(N,1),V(1:P+1,:)));
    int=(Kron.*kron(ones(N,1),V(1:P+1,:)))*ones(length(qW),1);
end