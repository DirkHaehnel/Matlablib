function [z,f] = hongfun(p, t, c, n);

c = c(:)';
kon = p(1);
koff = p(2);
kcl = p(3);
coef = ones(1,n);
for j=1:length(c)
    z0 = zeros(2*n, 1);
    z0(1) = 1;
    M = [-c(j)*kon*eye(n) koff*eye(n); c(j)*kon*eye(n) (-(koff+kcl)*eye(n) + kcl*diag(ones(1,n-1),-1))];
    for k=1:length(t)
        z(:,k) = expm(t(k)*M)*z0;
        f(k,j) = [coef coef]*z(:,k);
    end;
end
plot(t,f); drawnow

