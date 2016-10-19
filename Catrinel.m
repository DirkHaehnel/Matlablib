% Defoc Imaging Catrinel

close all

rhov = [0 5];
z = 0;

load metals

%nsi = 3.8732 + i*0.014529; % refractive index of Si
nsi = gold(wavelength==670);
nsio = 1.46; % refractive index of SiO
nair = 1; % refractive index of air
dsio = 10e-3;

focpos = 0.5;
mag = 100;
NA = 1.4;
lambda = 0.67;
n1 = 1;
n = 1;
n2 = [nsio nsi];
d1 = [];
d = 0;
d2 = dsio;

[intx, inty, intz, rho, phi, fxx0, fxx2, fxz, byx0, byx2, byz] = SEPDipole(rhov, z, NA, n1, n, n2, d1, d, d2, lambda, mag, focpos);
cnt=1;
for theta=0:10:90
    subplot(2,5,cnt);
    int = SEPImage(theta/180*pi,0,30,10,rho,fxx0,fxx2,fxz,byx0,byx2,byz);
    mim(int)
    title(['\theta = ' int2str(theta) '°'],'fontsize',12);
    cnt = cnt + 1;
end

