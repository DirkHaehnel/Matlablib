% DefocHeight

load metals;

focpos = 1;
rhov = [0 2];
z = 0.07;
NA = 1.35; 
n1 = 1.51;
n = 1.49;
% n2 = aluminum(wavelength==605);
n2 = 1;
d1 = [];
d = 0.1;
d2 = [];
lambda = 0.605;
mag = 200;

[intx, inty, intz, rho, phi, fxx0, fxx2, fxz, byx0, byx2, byz] = ...
    SEPDipole(rhov, z, NA, n1, n, n2, d1, d, d2, lambda, mag, focpos);
cnt=1;
for theta=0:10:90
    subplot(2,5,cnt);
    int = SEPImage(theta/180*pi,0,30,10,rho,fxx0,fxx2,fxz,byx0,byx2,byz);
    mim(int)
    title(['\theta = ' int2str(theta) '°'],'fontsize',12);
    cnt = cnt + 1;
end

