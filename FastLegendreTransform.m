function Y = FastLegendreTransform(lp,T)
    %p == the degree of the Legendre polynomial 
    %T == the domain, a=T(1)<...<T(end)=b
    %Returns Y, which satisfies Y(i)=Lp(ksi(T(i))), where Lp is the 
    %legendre polynomial of degree p,and ksi is the transformation 
    %ksi(t)=(2*t-(a+b))/(b-a), -1<=ksi<=1.
    %
    %The inverse transformation is t(ksi)=a*(1-ksi)/2+b*(1+ksi)/2
    a=T(1);
    b=T(end);
    Y=lp((2.*T-(a+b))./(b-a));
    
end