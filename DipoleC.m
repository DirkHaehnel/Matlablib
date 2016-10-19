function [fr,ff,fz] = DipoleC(theta,nmax,drho,n0,n,n1,d0,d,d1)

% [fr,ff,fz] = DipoleC(theta,nmax,drho,n0,n,n1,d0,d,d1) calculates the electric 
% field amplitudes of dipole radiation along emission angle theta 
% of a dipole at distance drho from an interface within a layered cylindrical structure 
% theta - polar angle of radiation 
% nmax - cut off of angular expansion
% drho  - molecule's distance from the inner side of its cylindrical layer
% n0 - vector of refractive indices of the inner stack of cylindrical layers
% n  - refracive index of the molecule's cylindrical layer
% n1 - vector of refractive indices of the outer stack of cylindrical layers
% d0 - vector of layer thickness values of the inner stack of cylindrical layers ( length(d0)=length(n0)-1 )
% d  - thickness of molecule's cylindrical layer
% d1 - vector of layer thickness values of the outer stack of cylindrical layers ( length(d1)=length(n1)-1 )

bj = inline('besselj(n,x)','n','x');
bh = inline('besselj(n,x) + i*bessely(n,x)','n','x');
bjx = inline('0.5*(besselj(n-1,x) - besselj(n+1,x))','n','x');
bhx = inline('0.5*(besselj(n-1,x) - besselj(n+1,x) + i*bessely(n-1,x) - i*bessely(n+1,x))','n','x');

n0 = n0(:).'; 
n1 = n1(:).';
d0 = d0(:).'; 
d1 = d1(:).';
nv = [n0 n n1];
rv = cumsum([d0 d d1]);
w = nv(end)*cos(theta(:));
q = sqrt(n^2 - w.^2);
r = sum(d0)+drho;

c0 = zeros(length(w),2*nmax+1,4*length(nv)-4);
fr = zeros(length(w),2*nmax+1,2);
ff = fr;
fz = fr;

% r-dipole
for m=-nmax:nmax
    c0(:,nmax+m+1,4*length(n1) - 1) = (-i*m./r.*bj(m,q.*r)./q.^2).';
    c0(:,nmax+m+1,4*length(n1)) = (-i*w.*bjx(m,q.*r)./q/n^2).';
    c0(:,nmax+m+1,4*length(n1) + 1) = (-i*m./r.*bh(m,q.*r)./q.^2).';
    c0(:,nmax+m+1,4*length(n1) + 2) = (-i*w.*bhx(m,q.*r)./q/n^2).';
end
fr = Cylindrical(w,-nmax:nmax,nv,rv,c0);

% phi-dipole
for m=-nmax:nmax
    c0(:,nmax+m+1,4*length(n1) - 1) = (-bjx(m,q.*r)./q).';
    c0(:,nmax+m+1,4*length(n1)) = (-m.*w./r.*bj(m,q.*r)./q.^2/n^2).';
    c0(:,nmax+m+1,4*length(n1) + 1) = (-bhx(m,q.*r)./q).';
    c0(:,nmax+m+1,4*length(n1) + 2) = (-m.*w./r.*bh(m,q.*r)./q.^2/n^2).';
end
ff = Cylindrical(w,-nmax:nmax,nv,rv,c0);

% z-dipole
for m=-nmax:nmax
    c0(:,nmax+m+1,4*length(n1) - 1) = 0;
    c0(:,nmax+m+1,4*length(n1)) = (bj(m,q.*r)/n^2).';
    c0(:,nmax+m+1,4*length(n1) + 1) = 0;
    c0(:,nmax+m+1,4*length(n1) + 2) = (bh(m,q.*r)/n^2).';
end
fz = Cylindrical(w,-nmax:nmax,nv,rv,c0);



