dist = 0.40;
rhofield = [0 3];
zfield = [10 20];
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
resolution = [40 10];
ring = [];
maxm = 10;

rmax = rhofield(end)/1.5;
[x,y] = meshgrid(-rmax:0.05:rmax,-rmax:0.05:rmax);
rr = sqrt(x.^2+y.^2);
pp = angle(x+i*y);

global pd
pd = 1/5e-5;

over = 1.5e3;
atf = '';
velo = 5e-5;
sat = 0;
tic
exc = DICExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focposexc, pow, atf, resolution, maxm);
mdf = DICExc2MDF(exc, NA, n0, n, n1, focposdet, lamem, mag, av, zpin, atf, kappa, lt, pulse, sat, triplet);
[modres, autotime] = FCS(exc.rho,exc.z,mdf.volx1,mdf.voly1,mdf.volx2,mdf.voly2,[],velo);

[w0, a0, vx] = simplex('GaussFcs',[0.4 0.1 1],[0 0 0],[],[],[],av/mag,[lamex lamem]/n,dist,autotime,modres(:,1:2),modres(:,3:4)/2,0);
dc = 1/pd/5e-5;

% save 2fFCSFlowXRes fullres autotime w0 a0 dc vx
toc
