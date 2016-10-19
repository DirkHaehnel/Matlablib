% PSFInterference

over = 5e3;
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
lamex = 0.64;
focpos = 0;
lamem = 0.67;
mag = 60;
av = 25;
zpin = 0e3;
atf = [];
kappa = 1;,
fd = 3e3;
resolution = [20 20];

exc = GaussExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, atf, resolution);

for shift=0:0.01:0.5
    exc1 = GaussExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos+shift, atf, resolution);
    exc1.fxc = exc1.fxc + exc.fxc;
    exc1.fxs = exc1.fxs + exc.fxs;
    exc1.fyc = exc1.fyc + exc.fyc;
    exc1.fys = exc1.fys + exc.fys;
    exc1.fzc = exc1.fzc + exc.fzc;
    exc1.fzs = exc1.fzs + exc.fzs;

    FocusImage2D(exc.rho,exc.z,cat(4,cat(3,exc1.fxc,exc1.fxs),cat(3,exc1.fyc,exc1.fys),cat(3,exc1.fzc,exc1.fzs)))
    text(0.2,0.1,['\delta = ' mint2str(shift*1e3,3) ' nm'],'units','normal','color','white')
    eval(['print -dpng tmp' mint2str(100*shift,2)])
end


return
mdf = GaussExc2MDF(exc, NA, n0, n, n1, focpos, lamem, mag, av, zpin, atf, kappa);
rho = exc.rho;
z = exc.z;
maxm = (size(mdf.volx,3)-1)/2;

