function [err, c, zz, z] = ExpPur(p, t, x);

t = t(:)-min(t); p = p(:)'; 
[m, n] = size(x);
if m<n x=x'; tmp=n; n=m; m=tmp; end
zz = exp(-t*p);
zz = zz./(ones(m,1)*sum(zz));

for j=1:n
    c(:,j) = lsqnonneg(zz,x(:,j));
    z(:,j) = zz*c(:,j);
end

semilogy(t, z, t, x, '.'); drawnow

err = sum(sum((x-z).^2./abs(z)))
%err = sum((x-z).^2)

