NA = 1.2;
fd = 3e3;
n0 = 1.333;
n = 1.333;
n1 = 1.333;
d0 = [];
d = [];
d1 = [];
lamex = 0.47;
over = 5e3;
focpos = 10;
av = 50/2; % confocal pinhole radius
lamem = 0.53;
mag = 60;
zpin = 0; 
atf = [];
resolution = 50;
kappa = 1; 
lt = [];
pulse = [];
sat = []; 
triplet = [];

rhofield = [0 1];
zfield = [8 12];

exc = GaussExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, atf, resolution);
wave = abs(exc.fxc(:,:,1)).^2 + abs(exc.fyc(:,:,1)).^2 + abs(exc.fzc(:,:,1)).^2 + sum(abs(exc.fxc(:,:,2:end)).^2 + abs(exc.fyc(:,:,2:end)).^2 + abs(exc.fzc(:,:,2:end)).^2,3)/2;
subplot(121)
mpcolor(exc.rho, exc.z, wave, 'lr');
xlabel('\rho [\mum]')
ylabel('\itz\rm [\mum]')
title('excitation intensity distribution');

mdf = GaussExc2MDF(exc, NA, n0, n, n1, focpos, lamem, mag, av, zpin, atf, kappa, lt, pulse, sat, triplet);
subplot(122)
FocusImage2D(exc.rho,exc.z,mdf.volx+mdf.voly);
xlabel('\rho [\mum]')
ylabel('\itz\rm [\mum]')
title('MDF');

figure

FocusImage3D(exc.rho,exc.z,mdf.volx+mdf.voly);

figure
a = [flipud(mdf.volx(:,213,1)+mdf.voly(:,213,1)); mdf.volx(:,213,1)+mdf.voly(:,213,1)]; 
b = mdf.volx(1,:,1)+mdf.voly(1,:,1); 
c = exc.rho(:,1)'*(mdf.volx(:,:,1)+mdf.voly(:,:,1)); 
a=a/max(a); b=b/max(b); c=c/max(c); 
plot([-flipud(exc.rho(:,1)); exc.rho(:,1)],a,exc.z(1,:)-10,b,exc.z(1,:)-10,c)
xlabel('distance (\mum');
ylabel('MDF')
grid
legend({'lateral distance','along optical axis','membrane scan'})

return

% calculation of depolarization by objective: 

kappa = 1;

% repeat above

plot(exc.z(1,:),sum(exc.rho.*mdf.volx(:,:,1)),exc.z(1,:),sum(exc.rho.*mdf.voly(:,:,1)))
max(sum(exc.rho.*mdf.volx(:,:,1))/sum(exc.rho.*mdf.voly(:,:,1)))

% normally: 3
theta = (0.5:1e3)/1e3*pi;
2*sum(sin(theta).*cos(theta).^4)/sum(sin(theta).*cos(theta).^2.*sin(theta).^2)

% in plane image
phi=(0:360)/360*2*pi;
fx=exc.fxc(:,1,1)*ones(size(phi));fy=exc.fyc(:,1,1)*ones(size(phi));fz=exc.fzc(:,1,1)*ones(size(phi));
for j=1:2 
    fx=fx+exc.fxc(:,1,j+1)*cos(j*phi)+exc.fxs(:,1,j)*sin(j*phi);
    fy=fy++exc.fyc(:,1,j+1)*cos(j*phi)+exc.fys(:,1,j)*sin(j*phi);
    fz=fz+exc.fzc(:,1,1+j)*cos(j*phi)+exc.fzs(:,1,j)*sin(j*phi);
end
pcolor(exc.rho(:,1)*cos(phi),exc.rho(:,1)*sin(phi),abs(fx).^2+abs(fy).^2+abs(fz).^2); axis image; shading flat


