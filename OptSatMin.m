function [err, z, c] = OptSatMin(p,x,y);

x = x(:);
y = y(:);
p = p(:)';

z = [ones(size(x)) x./polyval([1./p 1],abs(x))];

c = z\y;

z = z*c;

plot(x,y,'o',x,z); drawnow

err = sum((y-z).^2);