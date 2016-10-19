dist = 0.40;
rhofield = [0 1];
zfield = [12 18];
NA = 1.14;
fd = 3e3;
n0 = 1.333;
n = n0;
n1 = n0;
d0 = [];
d = 0;
d1 = [];
lamex = 0.64;
focposexc = [15 dist/2 0 0 0];
pow = 1;
lamem = 0.67;
mag = 60;
av = 100;
focposdet = 15;
zpin = 0e3;
kappa = 0;
lt = [];
pulse = [0.05 25]/2; % laser characteristics in units of lifetime
triplet = 0;
resolution = [30 10];
ring = [];
maxm = 10;

sat = 0;
over = 1.5e3;
atf = [];

close all

exc = DICExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focposexc, pow, atf, resolution, maxm);
mdf = DICExc2MDF(exc, NA, n0, n, n1, focposdet, lamem, mag, av);

FocusImage3D(exc.rho,exc.z,mdf.volx1+mdf.voly1)
FocusImage3D(exc.rho,exc.z,mdf.volx2+mdf.voly2,[],1)

