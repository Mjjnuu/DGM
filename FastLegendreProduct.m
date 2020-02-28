function Y = FastLegendreProduct(f,lp,x1,x2)
    F = @(w) f(x1.*(1-w)./2+x2.*(1+w)./2);
    %L = @(w) lp(x1.*(1-w)./2+x2.*(1+w)./2);
    %Y = L2norm(F,lp,-1,1);
    Y = FastL2norm(F,lp);
end    