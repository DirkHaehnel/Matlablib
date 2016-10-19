if ~exist('gold') load metals; end
NA = 1.4;
wd = 3e3;
tubelens = 180e3;
lamex = 0.635;
lamem = 0.67;
n0ex = [1.51 gold(wavelength==1e3*lamex)];
n0em = [1.51 gold(wavelength==1e3*lamem)];    
n = 1.333;
n1 = 1.333;
d = [];
d1 = [];
focpos = 0;
mag = tubelens/wd; %60;
av = 100/2; %pinhole radius!!!
zpin = 0; 
atf = [1.51 0];
kappa = 0; 
lt = 1;
sat = 0;
rhofield = [0 1];
zfield = [0 1];
over = [NA*wd 0 5e3];
ring = 0; %n*wd+1;

[volx0, voly0, rho, z] = MDF(rhofield, zfield, NA, n0ex([1 1]), n, n1, d0, d, d1, lamex, over, focpos, lamem, mag, av, zpin, atf, kappa, lt, sat, [], ring);

for j=1:10
    d0 = 0.005*j;
    [volx, voly, rho, z] = MDF(rhofield, zfield, NA, [n0ex n0em], n, n1, d0, d, d1, lamex, over, focpos, lamem, mag, av, zpin, atf, kappa, lt, sat, [], ring);

    pcolor([-flipud(rho);rho],[z;z],[flipud(volx0(:,:,1));volx(:,:,1)]); shading interp; axis image; set(gcf,'renderer','zbuffer')
    times16
    title(['film thickness = ' mnum2str(d0*1e3,3,0) ' nm']);
    eval(['print -dpng -r100 tmp' mint2str(j-1,1)]);
end