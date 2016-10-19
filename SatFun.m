function [err, c, z] = SatFun(p,x,y);

x = x(:);
y = y(:);

z = [ones(size(x)) x.^p(2)./(1+x.^p(2)/p(1))];

c = lsqnonneg(z,y);

z = z*c;

plot(x,y,x,z,'o'); drawnow

err = sum((y-z).^2);