function y=getPulse(N,x0,X,t,results)
    P=size(results, 1)/N-1;
    startT=t(1);
    endT=t(end);
    nT=length(t);
    deltaT=(endT-startT)/(nT-1);
    for e=1:N
        if X(e)<=x0 && x0<X(e+1)
            break
        end   
    end
    
    U=results((P+1)*(e-1)+1:(P+1)*e , :);

    ksi=(2*x0-(X(e+1)+X(e)))/(X(e+1)-X(e));
    LPfast=zeros(1,11);
    LPfast(1) = 1;
    LPfast(2) = ksi;
    LPfast(3) = 1/2*(3*ksi.^2-1);
    LPfast(4) = 1/2*(5*ksi.^3-3.*ksi);
    LPfast(5) = 1/8*(35*ksi.^4-30*ksi.^2+3);
    LPfast(6) = 1/8*(63*ksi.^5-70*ksi.^3+15*ksi);
    LPfast(7) = 1/16*(231*ksi.^6-315*ksi.^4+105*ksi.^2-5);
    LPfast(8) = 1/16*(429*ksi.^7-693*ksi.^5+315*ksi.^3-35*ksi);
    LPfast(9) = 1/128*(6435*ksi.^8-12012*ksi.^6+6930*ksi.^4-1260*ksi.^2+35);
    LPfast(10) = 1/128*(12155*ksi.^9-25740*ksi.^7+18018*ksi.^5-4620*ksi.^3+315*ksi);
    LPfast(11) = 1/256*(46189*ksi.^10-109395*ksi.^8+90090*ksi.^6-30030*ksi.^4+3465*ksi.^2-63);
    
    y=LPfast(1:P+1)*U(:,1:end-1);
    
end