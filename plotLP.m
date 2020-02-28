function  Y = plotLP(t,x,U)
    N=length(x)-1;
    P=size(U,1)-1;
    Y=zeros(1,length(t));
    
    LPfast=cell(11,1);
    LPfast{1} = @(x) 1;
    LPfast{2} = @(x) x;
    LPfast{3} = @(x) 1/2*(3*x.^2-1);
    LPfast{4} = @(x) 1/2*(5*x.^3-3.*x);
    LPfast{5} = @(x) 1/8*(35*x.^4-30*x.^2+3);
    LPfast{6} = @(x) 1/8*(63*x.^5-70*x.^3+15*x);
    LPfast{7} = @(x) 1/16*(231*x.^6-315*x.^4+105*x.^2-5);
    LPfast{8} = @(x) 1/16*(429*x.^7-693*x.^5+315*x.^3-35*x);
    LPfast{9} = @(x) 1/128*(6435*x.^8-12012*x.^6+6930*x.^4-1260*x.^2+35);
    LPfast{10} = @(x) 1/128*(12155*x.^9-25740*x.^7+18018*x.^5-4620*x.^3+315*x);
    LPfast{11} = @(x) 1/256*(46189*x.^10-109395*x.^8+90090*x.^6-30030*x.^4+3465*x.^2-63);
    
    for element = 1:N
        [subdomain,l,u] = GetDomain2(t,x(element),x(element+1));
        values=zeros(1,length(subdomain));
        for p=0:P
            values=values + FastLegendreTransform(LPfast{p+1},subdomain)*U(p+1,element);
        end    
        Y(l:u)=values;
    end    

end