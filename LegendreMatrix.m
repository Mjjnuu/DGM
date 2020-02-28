function LP=LegendreMatrix(qP)
    %Returns a 11xlength(qP) matrix containing the values of legendre
    %poynomials of degree 0-10 evaluated at quadrature points qP
    LP=zeros(11,length(qP));
    LP(1,:)=1;
    LP(2,:)=qP;
    LP(3,:)=1/2*(3*qP.^2-1);
    LP(4,:)=1/2*(5*qP.^3-3.*qP);
    LP(5,:)=1/8*(35*qP.^4-30*qP.^2+3);
    LP(6,:)=1/8*(63*qP.^5-70*qP.^3+15*qP);
    LP(7,:)=1/16*(231*qP.^6-315*qP.^4 ...
        +105*qP.^2-5);
    
    LP(8,:)=1/16*(429*qP.^7-693*qP.^5 ...
        +315*qP.^3-35*qP);
    
    LP(9,:)=1/128*(6435*qP.^8-12012*qP.^6 ...
        +6930*qP.^4-1260*qP.^2+35);
    
    LP(10,:)=1/128*(12155*qP.^9-25740*qP.^7 ...
        +18018*qP.^5-4620*qP.^3+315*qP);
    
    LP(11,:)=1/256*(46189*qP.^10-109395*qP.^8 ...
        +90090*qP.^6-30030*qP.^4+3465*qP.^2-63);
end