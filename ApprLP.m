function U = ApprLP(f,x,P) %ApprLP(f,t,x,P)
    N=length(x)-1;
    
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

    diagonal=2./(2*(0:P)+1);
    Diagonal=zeros(1,(P+1)*N);

    %(P+1)xN matrix containing integrals F_{ij}=\int_j f*L_i dx
    F=zeros(P+1,N);

    for i=1:N
        for power = 0:P
            F(power + 1, i)=FastLegendreProduct(f,LPfast{power+1},x(i),x(i+1));
        end 
        Diagonal((i-1)*(P+1)+1  : i*(P+1))=diagonal;
    end

    F=reshape(F,[],1);

    %Mass matrix
    I=(1:(P+1)*N);
    M=sparse(I, I, Diagonal);
    U=M\F;
    U=reshape(U,P+1,N);

end