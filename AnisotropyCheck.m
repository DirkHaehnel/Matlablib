% AnisotropyCheck

if 0
    [x,x,c]=pt3v2read('Y:\Anastasia\2009-01-13 Anisotropy\Laser1+2.pt3',[1 inf]);
    y3a=mhist(x(c==1),t3);
    y3b=mhist(x(c==3),t3);
    [x,x,c]=pt3v2read('Y:\Anastasia\2009-01-13 Anisotropy\Laser2.pt3',[1 inf]);
    y2a=mhist(x(c==1),t2);
    y2b=mhist(x(c==3),t2);
    [x,x,c]=pt3v2read('Y:\Anastasia\2009-01-13 Anisotropy\Laser1.pt3',[1 inf]);
    y1a=mhist(x(c==1),t1);
    y1b=mhist(x(c==3),t1);
    m1a=max(y1a);
    m1b=max(y1b);
    m2b=max(y2b);
    m2a=max(y2a);
else
    m1a = 63391;
    m1b = 24712;
    m2b = 57417;
    m2a = 24584;
end

aniso_ex = sqrt(m1a*m2b/m1b/m2a)

NA = 1.2;
fd = 3e3;
n0 = 1.333;
n = 1.333;
n1 = 1.333;
d0 = [];
d = [];
d1 = [];
lamex = 0.64;
over = [0 1.6e3];
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

rhofield = [0 5];
zfield = [0 20];

exc = GaussExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, atf, resolution);
mdf = GaussExc2MDF(exc, NA, n0, n, n1, focpos, lamem, mag, av, zpin, atf, kappa, lt, pulse, sat, triplet);

aniso_th = sum(sum(exc.rho.*mdf.volx(:,:,1)))/sum(sum(exc.rho.*mdf.voly(:,:,1)))

ind = exc.rho(:,1)<0.5; 
p = Simplex('Gauss',[0.2 0],[],[],[],[],[-flipud(exc.rho(ind,1)); exc.rho(ind,1)],[flipud(mdf.volx(ind,round(end/2),1)); mdf.volx(ind,round(end/2),1)],[],[],1);

2*p(2)
