%clear all
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
focpos = 10;
lamem = 0.47;
mag = 60;
av = 100;
zpin = 0e3;
atf = [];
kappa = 0;
lt = [];
fd = 3e3;
resolution = [30 10];

NatConst
extinction = 80e3;
tau = 4;
pulse = [1 25]/tau;
power = (5:10:3e2)*1e-6;
wavelength = 400e-9;
sigma = 1e3*extinction*log(10)/AvogadroConstant*1e8; % mum^2
satv = sigma*power/(PlanckConstant*SpeedOfLight/wavelength)*tau*1e-9;

rhofield = [0 1];
zfield = [8 12];

exc = GaussExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, atf, resolution);
for j=1:length(satv)
    mdf(j) = GaussExc2MDF(exc, NA, n0, n, n1, focpos, lamem, mag, av, zpin, atf, kappa, lt, pulse, satv(j));
    intx(j) = sum(exc.rho(:,1)'*mdf(j).volx(:,:,1));
    inty(j) = sum(exc.rho(:,1)'*mdf(j).voly(:,:,1));    
end

plot(1e6*power,1e9*intx,1e6*power,1e9*inty);
xlabel('excitation power (\muW)')
ylabel('fluorescence intensity (a.u.)')

figure
subplot(121); pcolor(exc.rho,exc.z,mdf(1).volx(:,:,1)); axis image; shading interp; colorbar; subplot(122); pcolor(exc.rho,exc.z,mdf(end).volx(:,:,1)); axis image; shading interp; colorbar; colormap hot