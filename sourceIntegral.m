function int = sourceIntegral(U,S)
    %U is vector of coefficients in Legendre basis in one element
    %S is the source function
    %Calculates integrals \int_I S(u_h)L_q dx, where I=[-1,1], 
    %u_h=\sum_{p=0}^P U(p)L_p, and L_i is the Legendre polynomial
    %of degree i \in {0,1,...,P}.
    P=length(U)-1;
    
    %Gauss-Legendre quadrature points and weights 
    quadWeights = [0.0812743883615744, 0.1806481606948574, 0.2606106964029354, ...
        0.3123470770400029, 0.3302393550012598, 0.3123470770400029, ...
        0.2606106964029354, 0.1806481606948574, 0.0812743883615744];

%     quadWeights = [0.0812743883615744; 0.1806481606948574; 0.2606106964029354; ...
%         0.3123470770400029; 0.3302393550012598; 0.3123470770400029; ...
%         0.2606106964029354; 0.1806481606948574; 0.0812743883615744];
    
    quadPoints = [-0.9681602395076261, -0.8360311073266358, -0.6133714327005904, ...
        -0.3242534234038089, 0,0.3242534234038089, ...
        0.6133714327005904, 0.8360311073266358, 0.9681602395076261];
    
    %4th order Gauss-Legendre
%     qW=[0.3478548451374538; 0.6521451548625461; 0.6521451548625461; 0.3478548451374538];
%     qP=[-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526];
%     LP=zeros(11,4);
%     LP(1,:)=1;
%     LP(2,:)=qP;
%     LP(3,:)=1/2*(3*qP.^2-1);
%     LP(4,:)=1/2*(5*qP.^3-3.*qP);
%     LP(5,:)=1/8*(35*qP.^4-30*qP.^2+3);
%     LP(6,:)=1/8*(63*qP.^5-70*qP.^3+15*qP);
%     LP(7,:)=1/16*(231*qP.^6-315*qP.^4 ...
%         +105*qP.^2-5);
%     
%     LP(8,:)=1/16*(429*qP.^7-693*qP.^5 ...
%         +315*qP.^3-35*qP);
%    
%     LP(9,:)=1/128*(6435*qP.^8-12012*qP.^6 ...
%         +6930*qP.^4-1260*qP.^2+35);
%     
%     LP(10,:)=1/128*(12155*qP.^9-25740*qP.^7 ...
%         +18018*qP.^5-4620*qP.^3+315*qP);
%    
%     LP(11,:)=1/256*(46189*qP.^10-109395*qP.^8 ...
%         +90090*qP.^6-30030*qP.^4+3465*qP.^2-63);
    
    
    
    
    
    
    %Calculating the values of Legendre polynomials of degree 0 to 10 at
    %9th degree Gauss-Legendre points
    LPvalues=zeros(11,9);
    LPvalues(1,:)=1;
    LPvalues(2,:)=quadPoints;
    LPvalues(3,:)=1/2*(3*quadPoints.^2-1);
    LPvalues(4,:)=1/2*(5*quadPoints.^3-3.*quadPoints);
    LPvalues(5,:)=1/8*(35*quadPoints.^4-30*quadPoints.^2+3);
    LPvalues(6,:)=1/8*(63*quadPoints.^5-70*quadPoints.^3+15*quadPoints);
    LPvalues(7,:)=1/16*(231*quadPoints.^6-315*quadPoints.^4 ...
        +105*quadPoints.^2-5);
    
    LPvalues(8,:)=1/16*(429*quadPoints.^7-693*quadPoints.^5 ...
        +315*quadPoints.^3-35*quadPoints);
   
    LPvalues(9,:)=1/128*(6435*quadPoints.^8-12012*quadPoints.^6 ...
        +6930*quadPoints.^4-1260*quadPoints.^2+35);
    
    LPvalues(10,:)=1/128*(12155*quadPoints.^9-25740*quadPoints.^7 ...
        +18018*quadPoints.^5-4620*quadPoints.^3+315*quadPoints);
   
    LPvalues(11,:)=1/256*(46189*quadPoints.^10-109395*quadPoints.^8 ...
        +90090*quadPoints.^6-30030*quadPoints.^4+3465*quadPoints.^2-63);
    
    %9th degree Gauss-Legendre quadrature
    %\sum_{i=1}^9 w_i * S(u_h(x_i))*L_q(x_i),
    %where w_i are the quadWeights and x_i quadrature points.
    %u_h=\sum_{p=0}^P U(p)L_p
    for q=0:P
        int(q+1,1)=sum(quadWeights.*(S( U'*LPvalues(1:P+1,:) ) .* LPvalues(q+1,:)) );
        %int(q+1,1)=sum((S( U'*LPvalues(1:P+1,:) ) .* LPvalues(q+1,:))*quadWeights );
        
        %int(q+1,1)=sum((S( U'*LP(1:P+1,:) ) .* LP(q+1,:))*qW );
        %int(q+1,1)=(S( U'*LP(1:P+1,:) ) .* LP(q+1,:))*qW ;
    end    
    %int(1:P+1,1)=(S( U'*LP(1:P+1,:) ) .* LP(1:P+1,:))*qW;
    
    
    %int=(S( U'*LP ) .* LP)*qW;
    
        
end