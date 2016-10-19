function [err, c, z] = ExpPolyFun(p, t, x)

t = t(:); x = x(:); p = 1./p(:)';
for j=1:length(p)
    z(:,j) = exp(-t.^j*p(j));
    z(:,j) = z(:,j)./sum(z(:,j));
end

c = lsqnonneg([ones(length(x),1) z],x);
z = [ones(length(x),1) z]*c;

plot(t, z, t, x, 'o'); drawnow

err = sum((x-z).^2./abs(z))

