function [err, c, z, A] = ExpGauss(p,t,y);

if isempty(t)
	t = 1:length(y);
end
t = t(:);
y = y(:);

A = [exp(-t*p(1)) exp(-(t-p(2)).^2/p(3)^2/2)];
c = lsqnonneg(A,y);
z = A*c;
plot(t,y,'o',t,z); drawnow
err = sum((z-y).^2./abs(z));

