% program PhaseHalfPlateExc

rhofield = [0 1];
zfield = [-2 2];
NA = 1.2;
fd = 3e3;
n0 = 1.33;
n = 1.33;
n1 = 1.33;
d0 = [];
d = 0;
d1 = [];
lamex = 0.514;
over = 5e3;
focpos = 0;
atf = [];
resolution = [25 5];
maxm = 50;

ring = '-i*pi.*(psi<pi)';
[fxc, fxs, fyc, fys, fzc, fzs, rho, z] = GaussExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, atf, resolution, ring, maxm);
FocusImage3D(rho,z,cat(4,cat(3,fxc,fxs),cat(3,fyc,fys),cat(3,fzc,fzs)))
