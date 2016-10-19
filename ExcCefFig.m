rhofield=[0 3]; zfield=[0 10]; NA=1.4; n0=1.51; n=1.333; n1=1.333; d0=[]; d=[]; d1=[]; 
lamex=0.514; lamem = 0.57;
over=[3e3*1.4 2e6 5e3]; atf=[];
mag = 100; av = 100/2; iso = []; zpin = 0;
fpv = 0:0.2:10;

for j=1:length(fpv)
    focpos = fpv(j);
    [exc, excx, excz, rho, z] = NewExc(rhofield, zfield, NA, n0, n, n1, d0, d, d1, lamex, over, focpos, atf);
    res = NewCEF(rho, z, NA, n0, n, n1, d0, d, d1, lamem, mag, focpos, av, zpin, atf, iso);
    pcolor([-flipud(rho); rho], [z; z], [flipud(res/max(res(:))); exc/max(exc(:))]); 
    shading interp; axis image; times16;
    set(gca,'xticklabel',int2str(abs(str2num(get(gca,'xticklabel')))));
    xlabel('\rho [\mum]');
    ylabel('\itz \rm[\mum]');
    eval(['print -djpeg bld' mint2str(j,2)])
end

