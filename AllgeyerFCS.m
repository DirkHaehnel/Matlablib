NA = 1.2;
fd = 3e3;
n0 = 1.333;
n = 1.333;
n1 = 1.333;
d0 = [];
d = [];
d1 = [];
lamex = 0.47;
over = 3e3;
focpos = 10;
av = 50/2; % confocal pinhole radius
lamem = 0.53;
mag = 60;
zpin = 0; 
atf = [];
resolution = 50;
kappa = []; 
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
a=a/max(a); b=b/max(b); c=c/max(c); 
plot([-flipud(exc.rho(:,1)); exc.rho(:,1)],a,exc.z(1,:)-10,b)
xlabel('distance (\mum');
ylabel('MDF')
grid
legend({'lateral distance','along optical axis'})

figure 

[modres, autotime] = FCS(exc.rho,exc.z,mdf.volx,mdf.voly);
p = Simplex('Rigler', [0.001 0.01],[0 0],[],[],[], autotime, modres, [], 1);
for j=1:10 
    p = Simplex('Rigler', p,[0 0],[],[],[], autotime, modres, [], 1); 
end
axis([1 1e6 0 1.01])
xlabel('time (a.u.)'); ylabel('normalized ACF')
text(1e4,0.8,{'\tau_1 = 227','\tau_2 = 8150'})



