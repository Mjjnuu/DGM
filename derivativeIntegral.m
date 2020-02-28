function int = derivativeIntegral(N,U,F,qW,qP,LPminus,LPplus,V,epsilon)
    %INPUTS:
    
    %U is vector of coefficients in Legendre basis 
    %U=[U_1;U_2;...;U_N], where U_i is the coefficient vector on element i.
    
    %F(U(x)) is a function
    
    %qW is a row vector containing quadrature weigths
    %qP is a row vector containing quadrature points
    %|qW|=|qP|=q
    
    %LPminus is a (P+1) x q matrix. LPminus(p,i) is the value of the Legendre 
    %polynomial of degree p-1 evaluated at i:th quadrature point -epsilon i.e. 
    %LP(p,i)=L_(p-1)(qP(i)-epsilon), where L_i is the
    %Legendre polynomial of degree i-1.
    
    %V is a (P+1) x q matrix. v_i , i \in {0,1,...,P} are the testing
    %functions, degree(v_i)=i. V(p+1,i) is the value of p+1:th testing
    %function evaluated at i:th quadrature point, i.e. V(p+1,i)=v_p(qP(i)).
    
  
    %PURPOSE:
    %Calculates the integrals \int_I dF(u_h(x))/dx*v_j(x) dx using numerical 
    %quadrature determined by qW and qP. Here u_h=\sum_{p=0}^P U_ep*L_p on 
    %element e, U_e is the coefficient vector with respect to Legendre 
    %polynomial basis on element e.  
    
    %int is a N*(P+1)x1 matrix. int=[I_1;I_2;...;I_N], where I_e are the
    %integral vectors of element e. I_e(j)=\int_e dF(u_h(x))/dx*v_j(x) dx =
    %\sum_{i=1}^q qW(i)*F(u_h(qP(i)))*v_j(qP(i))
    P=length(U)/N-1;
    symmetricDifference=F(reshape(U,P+1,N)'*LPplus(1:P+1,:))-F(reshape(U,P+1,N)'*LPminus(1:P+1,:));
    symmetricDifference=symmetricDifference/(2*epsilon);
    
    int=(kron(symmetricDifference*diag(qW),ones(P+1,1)).*kron(ones(N,1),V(1:P+1,:)))*ones(length(qP),1);
end    
    
    
    