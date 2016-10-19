% DICPSFTest

rhofield = [0 3];
zfield = [-4 4];
NA = 1.14;
n0 = 1.333;
n = n0;%1;
n1 = n1;%1;
d0 = [];
d = 0;
d1 = [];
lamex = 0.64;

over = [0 1500 1];
dfoc = 0.2;
focpos = [[0 dfoc 0 0 0];[0 -dfoc 0 0 0]];
pow = 1;
lamem = 0.67;
mag = 60;
av = 100;
zpin = 0e3;
atf = [1.52 0];
kappa = 1;,
lt = [];
sat = 0;
fd = 3e3;
resolution = [20 lamex/0.2];
ring = [];
maxm = 10;

[volx1, voly1, volx2, voly2, rho, z, fxc1, fxs1, fyc1, fys1, fzc1, fzs1, fxc2, fxs2, fyc2, fys2, fzc2, fzs2, qv, qp] = DICPSF(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, pow, lamem, mag, av, zpin, atf, kappa, lt, sat, resolution, ring, maxm);
[al,ar] = FocusImage2D(rho,z,cat(3,fxc1,fxs1)); [bl,br] = FocusImage2D(rho,z,volx1);
a1 = abs([flipud(al);ar]).^2; b1 = [flipud(bl);br]; r = [-rho(end:-1:1,1);rho(:,1)];
kk = 10; plot(r,a1(:,kk)/max(a1(:,kk)),r,b1(:,kk)/max(b1(:,kk))); grid

return

[volx, voly, rho, z, fxc, fxs, fyc, fys, fzc, fzs, qv, qp] = FullMDF(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos(1,:), lamem, mag, av, zpin, atf, kappa, lt, sat, resolution, ring, maxm);
[al,ar] = FocusImage2D(rho,z,cat(3,fxc,fxs)); [bl,br] = FocusImage2D(rho,z,volx);
a = abs([flipud(al);ar]).^2; b = [flipud(bl);br]; 
kk = 1; plot(r,a(:,kk)/max(a(:,kk)),r,b(:,kk)/max(b(:,kk))); grid

return 

kk = 1; plot(r,a1(:,kk)/max(a1(:,kk)),r,b1(:,kk)/max(b1(:,kk))); grid
kk = 1; plot(r,a1(:,kk),r,a(:,kk)); grid
kk = 1; plot(r,b1(:,kk),r,b(:,kk)); grid




