function [err, c, zz, z] = expfunfun(p, t, x, y);

t = t(:)-min(t); x = x(:); p = p(:)';
zz = exp(-t*p);
zz = zz./(ones(length(x),1)*sum(zz));
col = ones(length(x),1)/length(x);

c = lsqnonneg([col zz y],x);
z = [col zz y]*c;

semilogy(t, z, t, x, '.'); drawnow

err = sum((x-z).^2./abs(z))
%err = sum((x-z).^2)

