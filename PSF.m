% PSF

over = 1.4e3;
cov = 0;
rhofield = [0 1];
zfield = [-3 3];
NA = 1.2;
n0 = 1.333;
n = 1.333;
n1 = 1.333;
d0 = [];
d = 0;
d1 = [];
lamex = 0.47;
focpos = 0;
lamem = 0.55;
mag = 60;
av = 100;
zpin = 0e3;
atf = [];
kappa = 1;
lt = 1;
fd = 3e3;
resolution = [30 10];

% load metals
% dau = 48/1e3;

% exc = GaussExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, atf, resolution);
% mdf = GaussExc2MDF(exc, NA, n0, n, n1, focpos, lamem, mag, av, zpin, atf, kappa, lt);
% nau = gold(wavelength==1e3*lamex);
% exc = GaussExc(rhofield, zfield, NA, fd, [n0 nau], n, n1, dau, d, d1, lamex, over, focpos, atf, resolution);
% nau = gold(wavelength==1e3*lamem);
% mdf = GaussExc2MDF(exc, NA, [n0 nau], n, n1, focpos, lamem, mag, av, zpin, atf, kappa, lt);

exc = GaussExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, atf, resolution);
mdf = GaussExc2MDF(exc, NA, n0, n, n1, focpos, lamem, mag, av, zpin, atf, kappa, lt);

rho = exc.rho(:,1);
z = exc.z;
maxm = size(exc.fxs,3);

figure
FocusImage2D(exc.rho,exc.z,mdf.volx+mdf.voly,[],'surf')

% phi = (0:pi/100:2*pi); row = ones(size(phi));
% ex = exc.fxc(:,1,1)*row; ey = exc.fyc(:,1,1)*row; ez = exc.fzc(:,1,1)*row;
% for j=1:maxm
%     ex = ex + exc.fxc(:,1,1+j)*cos(j*phi) + exc.fxs(:,1,j)*sin(j*phi);
%     ey = ey + exc.fyc(:,1,1+j)*cos(j*phi) + exc.fys(:,1,j)*sin(j*phi);
%     ez = ez + exc.fzc(:,1,1+j)*cos(j*phi) + exc.fzs(:,1,j)*sin(j*phi);    
% end
% mpcolor(rho,phi,abs(ex).^2)

