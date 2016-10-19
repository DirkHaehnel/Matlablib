function [err, c, z] = Poisson(m, x, y, flag)

m = m(:)';
x = x(:);
y = y(:);

if nargin>3 && ~isempty(flag)
    x = x/m(end);
    m = m(1:end-1);
end

if m<100
    z = exp(x*log(m))./(gamma(x+1)*exp(m));
else
    z = exp(x*log(m) - x.*log(x)+x - 0.5*log(2*pi*x) - m);
end
c = lsqnonneg(z,y);
plot(x,y,'o',x,z*c); drawnow
err = sum((y-z*c).^2);

