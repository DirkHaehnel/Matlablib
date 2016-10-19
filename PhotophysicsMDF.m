clear
close all

NA = 1.2;
fd = 3e3;
n0 = 1.333;
n = 1.333;
n1 = 1.333;
d0 = [];
d = [];
d1 = [];
lamex = 0.64;
over = 1.25e3;
focpos = 10;
av = 150/2; % confocal pinhole radius
lamem = 0.670;
mag = 60;
zpin = 0; 
atf = [];
resolution = 50;
kappa = 0; 
lt = [];
pulse = [];
sat = []; 
triplet = [];

rhofield = [0 1];
zfield = [10 14];

exc = GaussExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, atf, resolution);
wave = abs(exc.fxc(:,:,1)).^2 + abs(exc.fyc(:,:,1)).^2 + abs(exc.fzc(:,:,1)).^2;
mdf = GaussExc2MDF(exc, NA, n0, n, n1, focpos, lamem, mag, av, zpin, atf, kappa, lt, pulse, sat, triplet);

subplot(121);
FocusImage2D(exc.rho,exc.z,cat(4,cat(3,exc.fxc,exc.fxs),cat(3,exc.fyc,exc.fys),cat(3,exc.fzc,exc.fzs)),[],'vertical');
subplot(122);
FocusImage2D(exc.rho,exc.z,mdf.volx+mdf.voly,[],'vertical');

lt = [];
pulse = [];
sat = 0.1:0.1:2;

for j=1:length(sat)
    mdf = GaussExc2MDF(exc, NA, n0, n, n1, focpos, lamem, mag, av, zpin, atf, kappa, lt, pulse, sat, triplet);


return

figure
ef = squeeze(abs(exc.fxc(:,:,1)).^2+abs(exc.fyc(:,:,1)).^2+abs(exc.fzc(:,:,1)).^2); ef = ef/max(ef(:));
df = squeeze(mdf.volx(:,:,1) + mdf.voly(:,:,1)); df = df/max(df(:));
plot(exc.rho(:,1),ef(:,1),exc.rho(:,1),df(:,1))

figure
[b,a]=mHist(ef(:),0:0.02:1,df(:));
b = b/sum(b);
bar(a,b,1,'b')
disp([a*b sqrt(a.^2*b-(a*b)^2)])

figure
mim([fliplr([flipud(ef); ef]) [flipud(ef); ef]]',[fliplr([flipud(df); df]) [flipud(df); df]]')

