NA = 1.2;
fd = 3e3;
n0 = 1.333;
n = 1.333;
n1 = 1.333;
d0 = [];
d = [];
d1 = [];
lamex = 0.64;
over = [0 3.5e3];
focpos = 10;
av = 75/2; % confocal pinhole radius
lamem = 0.670;
mag = 60;
zpin = 0; 
atf = [];
resolution = 30;
kappa = 1; 
lt = [];
pulse = [];
sat = []; 
triplet = [];

rhofield = [0 2];
zfield = [5 15];

exc = GaussExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, atf, resolution);

tmp = exc;
tmp.fyc = 0*tmp.fyc;
tmp.fys = 0*tmp.fys;
tmp.fzc = 0*tmp.fzc;
tmp.fzs = 0*tmp.fzs;
mdf = GaussExc2MDF(tmp, NA, n0, n, n1, focpos, lamem, mag, av, zpin, atf, kappa, lt, pulse, sat, triplet);
fxx = sum(sum(exc.rho.*mdf.volx(:,:,1)));
fyx = sum(sum(exc.rho.*mdf.voly(:,:,1)));

tmp = exc;
tmp.fxc = 0*tmp.fxc;
tmp.fxs = 0*tmp.fxs;
tmp.fzc = 0*tmp.fzc;
tmp.fzs = 0*tmp.fzs;
mdf = GaussExc2MDF(tmp, NA, n0, n, n1, focpos, lamem, mag, av, zpin, atf, kappa, lt, pulse, sat, triplet);
fxy = sum(sum(exc.rho.*mdf.volx(:,:,1)));
fyy = sum(sum(exc.rho.*mdf.voly(:,:,1)));

tmp = exc;
tmp.fyc = 0*tmp.fyc;
tmp.fys = 0*tmp.fys;
tmp.fxc = 0*tmp.fxc;
tmp.fxs = 0*tmp.fxs;
mdf = GaussExc2MDF(tmp, NA, n0, n, n1, focpos, lamem, mag, av, zpin, atf, kappa, lt, pulse, sat, triplet);
fxz = sum(sum(exc.rho.*mdf.volx(:,:,1)));
fyz = sum(sum(exc.rho.*mdf.voly(:,:,1)));

