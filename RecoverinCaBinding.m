function x = RecoverinCaBinding(c,p,k1,k2)

for j=1:length(c)
    x(1,j) = fzero(@(z) RCBfun(z, c(j), p, k1, k2), [0 c(j)]);
    [err, y] = RCBfun(x(1,j),c(j),p,k1,k2);
    x(2,j) = y; 
end

