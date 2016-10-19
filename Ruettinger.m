% program Rüttinger

rhofield = [0 2];
zfield = [-3 3];
NA = 1.2;
fd = 3e3;
n0 = 1.333;
n = 1.333;
n1 = 1.333;
d0 = [];
d = [];
d1 = [];
lamex = 0.635;
over = [inf 2.37e3];
focpos = 20;
lamem = 0.67;
mag = 60;
av = 36.5/2; %pinhole radius!!!
zpin = 0;
atf = [];
kappa = 0;
lt = [];
satv = [0 0.1*exp((0:0.1:1)*log(100))];
resolution = [50 30];
ring = [];
maxm = 2;
pulse = [0.05 25]/2;
triplet = 0;

exc = GaussExc(rhofield, focpos + zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, atf, resolution, ring, maxm);
for j=1:length(satv)
    mdf = GaussExc2MDF(exc, NA, n0, n, n1, focpos, lamem, mag, av, zpin, atf, kappa, lt, pulse, satv(j), triplet);
    [modauto(:,j), modtime] = FCS(exc.rho,exc.z,mdf.volx,mdf.voly);
end

