function [err, c, z] = Exp2fun(p, t, x)

t = t(:); x = x(:); p = p(:)';
z = exp(-t.^2/2/p^2);

c = [ones(length(x),1) z]\x;
z = [ones(length(x),1) z]*c;

semilogx(t, z, t, x, 'o'); drawnow

err = sum((x-z).^2./abs(z))

