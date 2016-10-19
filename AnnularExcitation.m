% Annular Excitation

close all

rhofield = [0 2]; % radial region of calculation
zfield = 0; % in focus
NA = 1.3;
fd = 3e3; % principal plane focal distance
n0 = 1.51; % ref. ind. immersion medium
n = 1.33;
n1 = n;
d0 = [];
d = 0;
d1 = [];
lamex = 0.64; %excitation lambda
over = 5e3; % beam radius in mum
focpos = 0; % focus position with respect to surface
atf = [];
resolution = 50;
ring = '-1e3*(rad<0.9)';

[fxc, fxs, fyc, fys, fzc, fzs, rho, z] = GaussExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, atf, resolution, ring);

phi = 0:pi/200:2*pi;
ex = fxc(:,1,1)*cos(0*phi); ey = fyc(:,1,1)*cos(0*phi); ez = fzc(:,1,1)*cos(0*phi);
for j=1:2
    ex = ex + fxc(:,1,j+1)*cos(j*phi) + fxs(:,1,j)*sin(j*phi);
    ey = ey + fyc(:,1,j+1)*cos(j*phi) + fys(:,1,j)*sin(j*phi);
    ez = ez + fzc(:,1,j+1)*cos(j*phi) + fzs(:,1,j)*sin(j*phi);
end
ex = abs(ex).^2;
ey = abs(ey).^2;
ez = abs(ez).^2;
subplot(131)
pcolor(rho*cos(phi),rho*sin(phi),ex)
axis image; colormap hot; shading interp
subplot(132)
pcolor(rho*cos(phi),rho*sin(phi),ey)
axis image; colormap hot; shading interp
subplot(133)
pcolor(rho*cos(phi),rho*sin(phi),ez)
axis image; colormap hot; shading interp

