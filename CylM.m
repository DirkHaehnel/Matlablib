function [fr, ff] = CylM(m,q,r,ss,cc)

% [fr, ff] = CylM(m,q,r,ss) calculates the M-vector cylinder function
% of order m for the coordinates r with Bessel functions ss
% Output are the rho and phi-components, the z-component is always
% zero for M
% The output has to be completed by exp(i*(m*phi+w*z))

if nargin<4 || isempty(ss) || strcmp(ss,'j')
    bj = inline('besselj(n,x)','n','x');
    bjx = inline('0.5*(besselj(n-1,x) - besselj(n+1,x))','n','x');
elseif strcmp(ss,'h')
    bj = inline('besselj(n,x) + i*bessely(n,x)','n','x');
    bjx = inline('0.5*(besselj(n-1,x) - besselj(n+1,x) + i*bessely(n-1,x) - i*bessely(n+1,x))','n','x');
end

[mr, nr] = size(r);
[mq, nq] = size(q);
if mr==mq && nr==nq
    fr = i*m./r.*bj(m,q.*r);
    ff = -q.*bjx(m,q.*r);
elseif nr==1 && mq==1
    row = ones(mq,nq);
    col = ones(mr,nr);
    fr = i*m./(r*row).*bj(m,r*q);
    ff = -(col*q).*bjx(m,r*q);
elseif mr==1 && nq==1
    row = ones(mr,nr);
    col = ones(mq,nq);
    fr = i*m./(col*r).*bj(m,q*r);
    ff = -(q*row).*bjx(m,q*r);
else
    error('Wrong input');
end

if nargin>4 && ~isempty(cc)
    if strcmp(cc,'f')
        fr = ff;
    end
end