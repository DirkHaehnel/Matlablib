function y = AssociatedLegendrePoly(n,m,x,flag)

% Modeled after LegendrePoly.m by David Terr, Raytheon, 5-10-04

% Given nonnegative integer n and integer m, computes the 
% associated Legendre polynomial P_nm(x).
% if flag='s': P_nm(x)/sqrt(1-x^2)
% if flag='d': sqrt(1-x^2)*P_nm'(x)

m = abs(m);
if n==0
    y = ones(size(x)); 
else
    p = [1 0 -1];
    for j=2:n
        p = conv(p,[1 0 -1]);
    end
    for j=1:(n+m)
        p = (length(p)-(1:length(p))).*p;
        p(end) = [];
    end
    if nargin>3 & flag=='s'
        y = m*(-1)^m*sqrt(1-x.^2).^(m-1).*polyval(p/2^n/prod(1:n),x);
    elseif nargin>3 & flag=='d'
        y = -m*(-1)^m*x.*sqrt(1-x.^2).^(m-1).*polyval(p/2^n/prod(1:n),x);
        p = (length(p)-(1:length(p))).*p;
        p(end) = [];
        y = y + (-1)^m*sqrt(1-x.^2).^(m+1).*polyval(p/2^n/prod(1:n),x);
    else
        y = (-1)^m*sqrt(1-x.^2).^m.*polyval(p/2^n/prod(1:n),x);
    end
    y(isnan(y)) = 0;
end

