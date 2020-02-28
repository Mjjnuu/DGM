function norm = L2norm(a,b,x1,x2)
    %a and b are functions from [x1,x2] to real numbers.
    %Returns the L2-inner product of a and b.
    product = @(x) a(x).*b(x);
    
    norm = integral(product,x1,x2,'RelTol',0.001,'AbsTol',0.001);
    
end    