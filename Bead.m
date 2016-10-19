function [coef, M, V] = Bead(mv,nv,r,c0)

% coef = Bead(mv,nv,rv,c0) calculates the vector-spherical-harmonics
% coefficients coef for scattering the source field c0 at a spherical
% particle with radius r and refractive index nv(1)
% mv is a vector giving the principal order of the VSHs
%
% order of c0: from out-to-in: m,n/m,n 
% order of coef: from out-to-in: hm,hn/jm,jn 

bj = inline('besselj(n+0.5,x).*sqrt(pi./x/2)','n','x');
bh = inline('besselh(n+0.5,x).*sqrt(pi./x/2)','n','x');
bjx = inline('sqrt(pi/8./x)*(besselj(n+0.5,x)./x + besselj(n-0.5,x) - besselj(n+1.5,x))','n','x');
bhx = inline('sqrt(pi/8./x)*(besselh(n+0.5,x)./x + besselh(n-0.5,x) - besselh(n+1.5,x))','n','x');

coef = [];

n2 = nv(2); n1 = nv(1);
coef = [-((n2*bjx(mv,n2*r).*bj(mv,n1*r)-n1*bj(mv,n2*r).*bjx(mv,n1*r))./...
        (n2*bhx(mv,n2*r).*bj(mv,n1*r)-n1*bh(mv,n2*r).*bjx(mv,n1*r))).'.*c0(1,:); ...
        -((n1*bjx(mv,n2*r).*bj(mv,n1*r)-n2*bj(mv,n2*r).*bjx(mv,n1*r))./...
        (n1*bhx(mv,n2*r).*bj(mv,n1*r)-n2*bh(mv,n2*r).*bjx(mv,n1*r))).'.*c0(2,:)];

