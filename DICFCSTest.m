dist = 0.41;
rhofield = [0 1.5];
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
pulse = [0.05 12.5]/2;
triplet = 0;
resolution = [20 5];
ring = [];
maxm = 10;
sat = 0;
over = [0 1500 1];
atf = [];
aberration = [];

[modres, autotime, exc, mdf] = DICFCS(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focposexc, ...
    pow, lamem, mag, av, focposdet, zpin, atf, kappa, lt, pulse, sat, triplet, resolution, ring, aberration, maxm);

dif = 5e-5;

modres = 1./(4*dif*autotime + 0.5^2)./sqrt(4*dif*autotime + 5*0.5^2);
modres = [modres modres modres.*exp(-dist^2./(4*dif*autotime+0.5^2))];

close; tst=simplex('GaussFcs',[0.4 0.4],[0 0],[],[],[],1.2,100/60,1.33,0.67,autotime/1e5,modres(:,1:2),modres(:,3),0.41)
[err,pd] = GaussFcs(tst,1.2,100/60,1.33,0.67,autotime/1e5,modres(:,1:2),modres(:,3),0.41)