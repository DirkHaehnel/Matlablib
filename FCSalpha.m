function [alpha, alpha0] = FCSalpha(dif0, dif, c, a, b, t, w0, w);
% maximum likelihood estimator for FCS

dif = dif(:);
t = t(:)';
pp = sqrt(pi/2)^1.5*a^2*b;
alpha0 = c*FCSZfun(a,b,dif0,t) + c^2*FCSZfun(a,b,dif0).^2;
alpha = c*FCSZfun(a,b,dif,t) + c^2*FCSZfun(a,b,dif).^2;
if nargin>6
    alpha0 = alpha0/w0;
    alpha = alpha./(w*ones(1,length(t)));
    alpha = log((ones(length(dif),1)*alpha0)./alpha);
end


