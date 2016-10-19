clear all
NA = 1.2;
tubelens = 180e3;
n0 = 1.333;
n = 1.333;
n1 = 1.333;
d0 = [];
d = [];
d1 = [];
over = 1.7e3;
lamex = 0.405;
focpos = 0;
lamem = 0.47;
mag = 60;
av = 100;
zpin = 0e3;
atf = [];
kappa = 0;
lt = [];
fd = 3e3;
resolution = 20;

NatConst
extinction = 80e3;
tau = 4;
pulse = [1 25]/tau;
power = [10 20]*1e-6; %(5:10:3e2)*1e-6;
wavelength = 400e-9;
sigma = 1e3*extinction*log(10)/AvogadroConstant*1e8; % mum^2
satv = sigma*power/(PlanckConstant*SpeedOfLight/wavelength)*tau*1e-9;

rhofield = [0 1];
zfield = [-1-lamex/resolution/2 1];

exc = GaussExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, atf, resolution);
for j=1:length(satv)
    mdf(j) = GaussExc2MDF(exc, NA, n0, n, n1, focpos, lamem, mag, av, zpin, atf, kappa, lt, pulse, satv(j));
end

[fx,ex]= FocusImage2D(exc.rho,exc.z,mdf(1).volx/power(1));
im1 = [flipud(ex);fx];
[fx,ex]= FocusImage2D(exc.rho,exc.z,mdf(2).volx/power(2));
im2 = [flipud(ex);fx];
[fx,ex]= FocusImage2D(exc.rho,exc.z,mdf(1).volx/power(1)-mdf(2).volx/power(2));
im3 = [flipud(ex);fx];
CombineImages(cat(3,im1,im2,im3),3,1,'scale')


return % Fourier analysis

a = 100; b = 30;
subplot(121)
[e,f] = FocusImage2D(exc.rho,exc.z,mdf(1).volx);
w = fft2([flipud(e); f]);
w = abs(fftshift(w));
mim(w(end/2-a:end/2+a,end/2-b:end/2+b))

subplot(122)
[e,f] = FocusImage2D(exc.rho,exc.z,mdf(2).volx);
w = fft2([flipud(e); f]);
w = abs(fftshift(w));
mim(w(end/2-a:end/2+a,end/2-b:end/2+b))

colormap jet
