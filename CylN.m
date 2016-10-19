function [fr, ff, fz] = CylN(m,q,w,r,ss,cc)

% [fr, ff, fz] = CylN(m,q,w,r,ss) calculates the N-vector cylinder function
% of order m for the coordinates r with Bessel functions ss
% Output are the rho, phi, and z-components
% The output has to be completed by exp(i*(m*phi+w*z))

if nargin<5 || isempty(ss) || strcmp(ss,'j')
    bj = inline('besselj(n,x)','n','x');
    bjx = inline('0.5*(besselj(n-1,x) - besselj(n+1,x))','n','x');
elseif strcmp(ss,'h')
    bj = inline('besselj(n,x) + i*bessely(n,x)','n','x');
    bjx = inline('0.5*(besselj(n-1,x) - besselj(n+1,x) + i*bessely(n-1,x) - i*bessely(n+1,x))','n','x');
end

[mr, nr] = size(r);
[mq, nq] = size(q);
if mr==mq && nr==nq
    ff = -m.*w./r.*bj(m,q.*r);
    fr = i*q.*w.*bjx(m,q.*r);
    fz = q.^2.*bj(m,q.*r);
elseif nr==1 && mq==1
    row = ones(mq,nq);
    col = ones(mr,nr);
    ff = -m.*(col*w)./(r*row).*bj(m,r*q);
    fr = i*(col*(q.*w)).*bjx(m,r*q);
    fz = (col*q.^2).*bj(m,r*q);
elseif mr==1 && nq==1
    row = ones(mr,nr);
    col = ones(mq,nq);
    ff = -m.*(w*row)./(col*r).*bj(m,q*r);
    fr = i*((q.*w)*row).*bjx(m,q*r);
	fz = (q.^2*row).*bj(m,q*r);
else
    error('Wrong input');
end

if nargin>5 && ~isempty(cc)
    if strcmp(cc,'f')
        fr = ff;
    elseif strcmp(cc,'z')
        fr = fz;
    end
end

