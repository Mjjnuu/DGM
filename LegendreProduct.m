function Y = LegendreProduct(f,p,x1,x2)
    F = @(w) f(x1.*(1-w)./2+x2.*(1+w)./2);
    L = @(w) legendreP(p,w);
    Y = L2norm(F,L,-1,1);
end    