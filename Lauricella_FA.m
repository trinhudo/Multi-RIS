function out = Lauricella_FA(a, bvec, cvec, xvec)

N = length(xvec);

syms t x
HypergeomProd = 1;

for n = 1:N
    HypergeomProd = HypergeomProd...
        * hypergeom(bvec(n),cvec(n),xvec(n)*t);
end

intFunction = exp(-t)*t^(a-1)*HypergeomProd;

% result(x) = vpaintegral(intFunction, t, 0, Inf)/gamma(x);
result(x) = vpaintegral(intFunction, t, [0,Inf], 'RelTol', 1e-32, 'AbsTol', 0)/gamma(x);

out = double(vpa(result(a)));
end