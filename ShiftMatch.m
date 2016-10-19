function [err, c, z] = ShiftMatch(p, x, y);

x = x(:);
n = length(x);
t = (1:n)';
if ~(size(y,1)==length(x))
    y = y';
end

for j=1:size(y,2)
    z(:,j) = (1-mod(p(j),1))*y(rem(rem(t-floor(p(j))-1, n)+n,n)+1,j) + ...
        mod(p(j),1)*y(rem(rem(t-ceil(p(j))-1, n)+n,n)+1,j);
end

z = [ones(n,1) z];
c = lsqnonneg(z,x);
z = z*c;    
semilogy(t,x,'o',t,z); drawnow
err = sum((x-z).^2);
