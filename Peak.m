function [err, c, z, A] = Peak(p,t,y);

y = y(:);
if isempty(t)
	t = 1:length(y);
end
if ~(length(t)==length(y))
    A = ones(length(y),2);
    A(:,2) = (y-p(1)).^p(2).*exp(-(y-p(1))/p(3));
    err = A*t;
else
    t = t(:);
    nn = length(p)/2;
    
    A = ones(length(t),2);
    A(:,2) = (t-p(1)).^p(2).*exp(-(t-p(1))/p(3));
    c = lsqnonneg(A,y);
    z = A*c;
    plot(t,y,'o',t,z); drawnow
    err = sum((z-y).^2./abs(z))
end