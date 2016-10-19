% program for calculating radially/axially polarized focus

close all

load metals

rhofield = [0 1];
NA = 1.2;
lamex = 0.514;
nglass = 1.518;
nag = silver(wavelength==1e3*lamex);
theta_max = asin(NA/nglass);
theta_min = 0; %60/180*pi;
fd = 5e3;
n0 = [nglass nag];
n = nglass;
n1 = [nag nglass];
d0 = 0.05;
d = 0.12;
d1 = 0.05;
zfield = [0 d];
over = 1e6;
focpos = 0;
atf = [];
resolution = 50;
maxm = 3;

% Fig. 4
ring = [];
[fxc, fxs, fyc, fys, fzc, fzs, rho, z] = GaussExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, atf, resolution, ring, maxm);
[fx, fx] = FocusImage2D(rho,z,cat(3,fxc,fxs));
[fy, fy] = FocusImage2D(rho,z,cat(3,fyc,fys));
[fz, fz] = FocusImage2D(rho,z,cat(3,fzc,fzs));
maxfield = sqrt(max(fx(:)+fy(:)+fz(:)));
subplot(1,3,1)
FocusImage2D(rho,z,cat(3,fxc,fxs)/maxfield)
colorbar('h')
set(gca,'fontsize',12)
axis square
subplot(1,3,2)
FocusImage2D(rho,z,cat(3,fzc,fzs)/maxfield)
colorbar('h')
set(gca,'fontsize',12)
axis square
subplot(1,3,3)
phi = 0:pi/100:2*pi;
row = ones(size(phi));
fx = fxc(:,end/2,1)*row; fy = fyc(:,end/2,1)*row; fz = fzc(:,end/2,1)*row;
for j=2:3
    fx = fx + fxc(:,end/2,j)*cos((j-1)*phi)+fxs(:,end/2,j-1)*sin((j-1)*phi);
    fy = fy + fyc(:,end/2,j)*cos((j-1)*phi)+fys(:,end/2,j-1)*sin((j-1)*phi);
    fz = fz + fzc(:,end/2,j)*cos((j-1)*phi)+fzs(:,end/2,j-1)*sin((j-1)*phi);    
end
ff = abs(fx).^2+abs(fy).^2+abs(fz).^2;
pcolor(rho(:,1)*cos(phi),rho(:,1)*sin(phi),ff/maxfield^2)
shading interp
colorbar('h')
set(gca,'fontsize',12)
axis square
colormap jet

return

print -dpng -r300 Fig4_Dmitri


% Fig.1 
ring = 'cos(psi)'; 
[fxc, fxs, fyc, fys, fzc, fzs, rho, z] = GaussExc(rhofield, [0 2], NA, fd, n, n, n, [], 2, [], lamex, over, 1, atf, resolution, ring, maxm);
[fxc2, fxs2, fyc2, fys2, fzc2, fzs2] = RotateEMField(fxc, fxs, fyc, fys, fzc, fzs, pi/2);
fxc = fxc + fxc2;
fxs = fxs + fxs2;
fyc = fyc + fyc2;
fys = fys + fys2;
fzc = fzc + fzc2;
fzs = fzs + fzs2;
figure
subplot(1,2,1)
FocusImage2D(rho,1-z,cat(3,fxc,fxs)/maxfield)
colorbar('h')
set(gca,'fontsize',12)
axis square
subplot(1,2,2)
FocusImage2D(rho,1-z,cat(3,fzc,fzs)/maxfield)
colorbar('h')
set(gca,'fontsize',12)
axis square
colormap jet
print -dpng -r300 Fig1_Dmitri

% Fig.3 
ring = 'cos(psi)'; 
[fxc, fxs, fyc, fys, fzc, fzs, rho, z] = GaussExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, atf, resolution, ring, maxm);
[fxc2, fxs2, fyc2, fys2, fzc2, fzs2] = RotateEMField(fxc, fxs, fyc, fys, fzc, fzs, pi/2);
fxc = fxc + fxc2;
fxs = fxs + fxs2;
fyc = fyc + fyc2;
fys = fys + fys2;
fzc = fzc + fzc2;
fzs = fzs + fzs2;
figure
subplot(1,3,1)
FocusImage2D(rho,z,cat(3,fxc,fxs)/maxfield)
colorbar('h')
set(gca,'fontsize',12)
axis square
subplot(1,3,2)
FocusImage2D(rho,z,cat(3,fzc,fzs)/maxfield)
colorbar('h')
set(gca,'fontsize',12)
axis square
ring = ['-sin(psi).*(rad>sin(' num2str(theta_min) '))']; maxm = 3;
[fxc, fxs, fyc, fys, fzc, fzs, rho, z] = GaussExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, atf, resolution, ring, maxm);
[fxc2, fxs2, fyc2, fys2, fzc2, fzs2] = RotateEMField(fxc, fxs, fyc, fys, fzc, fzs, -pi/2);
fxc = fxc + fxc2;
fxs = fxs + fxs2;
fyc = fyc + fyc2;
fys = fys + fys2;
fzc = fzc + fzc2;
fzs = fzs + fzs2;
subplot(1,3,3)
FocusImage2D(rho,z,cat(3,fyc,fys)/maxfield)
colorbar('h')
set(gca,'fontsize',12)
axis square
colormap jet
print -dpng -r300 Fig3_Dmitri

% Fig.5 
d = 0.14;
zfield = [0 d];
focpos = d/2;
ring = 'cos(psi)'; 
[fxc, fxs, fyc, fys, fzc, fzs, rho, z] = GaussExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, atf, resolution, ring, maxm);
[fxc2, fxs2, fyc2, fys2, fzc2, fzs2] = RotateEMField(fxc, fxs, fyc, fys, fzc, fzs, pi/2);
fxc = fxc + fxc2;
fxs = fxs + fxs2;
fyc = fyc + fyc2;
fys = fys + fys2;
fzc = fzc + fzc2;
fzs = fzs + fzs2;
figure
subplot(2,2,1)
FocusImage2D(rho,z,cat(3,fxc,fxs)/maxfield)
colorbar('h')
set(gca,'fontsize',12)
axis square
subplot(2,2,2)
FocusImage2D(rho,z,cat(3,fzc,fzs)/maxfield)
colorbar('h')
set(gca,'fontsize',12)
axis square
focpos = d/4;
[fxc, fxs, fyc, fys, fzc, fzs, rho, z] = GaussExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, atf, resolution, ring, maxm);
[fxc2, fxs2, fyc2, fys2, fzc2, fzs2] = RotateEMField(fxc, fxs, fyc, fys, fzc, fzs, pi/2);
fxc = fxc + fxc2;
fxs = fxs + fxs2;
fyc = fyc + fyc2;
fys = fys + fys2;
fzc = fzc + fzc2;
fzs = fzs + fzs2;
subplot(2,2,3)
FocusImage2D(rho,z,cat(3,fxc,fxs)/maxfield)
colorbar('h')
set(gca,'fontsize',12)
axis square
subplot(2,2,4)
FocusImage2D(rho,z,cat(3,fzc,fzs)/maxfield)
colorbar('h')
set(gca,'fontsize',12)
axis square
colormap jet
print -dpng -r300 Fig5_Dmitri
