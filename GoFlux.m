function y=GoFlux(F,a,b)
    %F(a,b)={min(a<=u<=b)F(u) , if a<=b
    %        max(b<=u<=a)F(u) , otherwise
    booleanMax = a>b;
    booleanMin = a<=b;
    y=booleanMax.*max(F(a+kron((0:10),(b-a)./10)), [], 2)+...
        booleanMin.*min(F(a+kron((0:10),(b-a)./10)), [], 2);
end
    