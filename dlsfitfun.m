function [err, c, y] = dlsfitfun(p, t, x);

t = t(:); x = x(:);
y = exp(-p(1)*t-p(2)*t.^2/2);

c = y\x;
y = c*y;

semilogx(t,y,t,x,'o'); drawnow;

err = sum((y-x).^2)

