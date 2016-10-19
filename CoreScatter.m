function coef = CoreScatter(V, k0, k1, k2, r1, r2);

% coef = CoreScatter(V, k0, k1, k2, r1, r2) calculates the coefficients of the spherical 
% harmonics for scattering of an electromagnetic field on a spherical core-shell particle;
% source field is given by matrix V
% k0, k1, k2 are wave numbers of exterior, shell, and core material
% r1, r2 are inner and outer radius of shell

% Order of the coefficients:
% [an(in) am(in) bn(sh) bm(sh) cn(sh) cm(sh) dn(ex) dm(ex)]

% Order of Matrix:
% [ENex = ENwall; ENwall = ENin; EMex = EMwall; EMwall = EMin;  
%  BNex = BNwall; BNwall = BNin; BMex = BMwall; BMwall = BMin]

coef = [];

bj = inline('besselj(n+0.5,x).*sqrt(pi./x/2)','n','x');
bh = inline('besselh(n+0.5,x).*sqrt(pi./x/2)','n','x');
bjx = inline('sqrt(pi/8./x)*(besselj(n+0.5,x)./x + besselj(n-0.5,x) - besselj(n+1.5,x))','n','x');
bhx = inline('sqrt(pi/8./x)*(besselh(n+0.5,x)./x + besselh(n-0.5,x) - besselh(n+1.5,x))','n','x');

for n=1:size(V,2)
    M = [bhx(n,k0*r2), 0, -bhx(n,k1*r2), 0, -bjx(n,k1*r2), 0, 0, 0; ...
            0, 0, bhx(n,k1*r1), 0, bjx(n,k1*r1), 0, -bjx(n,k2*r1), 0; ...
            0, -bh(n,k0*r2), 0, bh(n,k1*r2), 0, bj(n,k1*r2), 0, 0; ...
            0, 0, 0, -bh(n,k1*r1), 0, -bj(n,k1*r1), 0, bj(n,k2*r1); ...
            -k0*bh(n,k0*r2), 0, k1*bh(n,k1*r2), 0, k1*bj(n,k1*r2), 0, 0 ,0; ...
            0, 0, -k1*bh(n,k1*r1), 0, -k1*bj(n,k1*r1), 0, k2*bj(n,k2*r1), 0; ...
            0, k0*bhx(n,k0*r2), 0, -k1*bhx(n,k1*r2), 0, -k1*bjx(n,k1*r2), 0, 0; ...
            0, 0, 0, k1*bhx(n,k1*r1), 0, k1*bjx(n,k1*r1), 0, -k2*bjx(n,k2*r1) ...
        ];
    
    coef = [coef M\V(:,n)];
end

