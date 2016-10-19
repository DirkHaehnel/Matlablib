load c:\Joerg\Doc\Nanophotonics&Optics\Khoptyar\profile.mat

rhofield = [0 1];
zfield = 0;
NA = 1.25;
fd = 164.5e3/100;
n0 = 1.52;
n = 1.46;
n1 = 1;
d0 = [];
d = 0.05;
d1 = [];
lamex = 0.488;
over = 4e3;
focpos = 0;
atf = [];
resolution = 50;

eflag = true;
% eflag = false;

ring = ['cos(psi)']; maxm = 3;
[fxc, fxs, fyc, fys, fzc, fzs, rho, z] = GaussExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, atf, resolution, ring, maxm);
[fxc2, fxs2, fyc2, fys2, fzc2] = RotateEMField(fxc, fxs, fyc, fys, fzc, fzs, pi/2);

fxc = fxc + fxc2;
fzc = fzc + fzc2;
rho = [-flipud(rho); rho];
fx = [flipud(fx); fx];
fz = [flipud(fz); fz];

ff = [];
for phi=0:pi/4:pi/2
   ff = [ff squeeze(abs(sum(fxc(:,1,:)*sin(phi)^2+fzc(:,1,:)*cos(phi)^2,3)).^2)];
end
ff = [flipud(ff); ff];
ff = ff./(ones(size(ff,1),1)*max(ff));

rad = rad1; cnt = 1; r0v = []; err = [];
for r0 = 0.4:0.01:0.7
    fm = interp1(rho,ff,rad(:,1)-r0);
    c = lsqnonneg([ones(size(rad,1),1) fm],rad(:,2)-min(rad(:,2)));
    err(cnt) = sum((rad(:,2)-min(rad(:,2))-[ones(size(rad,1),1) fm]*c).^2);
    r0v(cnt) = r0;
    cnt = cnt+1;
end
r0 = r0v(err==min(err));
fm = interp1(rho,ff,rad(:,1)-r0);
c = lsqnonneg([ones(size(rad,1),1) fm],rad(:,2)-min(rad(:,2)));
plot(rad(:,1),rad(:,2)-min(rad(:,2)),rad(:,1),[ones(size(rad,1),1) fm]*c)
print -dpng -r300 RadialMode1

rad = rad2; cnt = 1; r0v = []; err = [];
for r0 = 0.4:0.01:0.7
    fm = interp1(rho,ff,rad(:,1)-r0);
    c = lsqnonneg([ones(size(rad,1),1) fm],rad(:,2)-min(rad(:,2)));
    err(cnt) = sum((rad(:,2)-min(rad(:,2))-[ones(size(rad,1),1) fm]*c).^2);
    r0v(cnt) = r0;
    cnt = cnt+1;
end
r0 = r0v(err==min(err));
fm = interp1(rho,ff,rad(:,1)-r0);
c = lsqnonneg([ones(size(rad,1),1) fm],rad(:,2)-min(rad(:,2)));
plot(rad(:,1),rad(:,2)-min(rad(:,2)),rad(:,1),[ones(size(rad,1),1) fm]*c)
print -dpng -r300 RadialMode2

ring = ['-sin(psi)']; maxm = 3;
[fxc, fxs, fyc, fys, fzc, fzs] = GaussExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, atf, resolution, ring, maxm);
[fxc2, fxs2, fyc2, fys2] = RotateEMField(fxc, fxs, fyc, fys, fzc, fzs, -pi/2);

fyc = fyc + fyc2;
ff = squeeze(abs(sum(fyc(:,1,:),3)).^2);
ff = [flipud(ff); ff]; ff = ff/max(ff);

rad = azi1; cnt = 1; r0v = []; err = [];
for r0 = 0.4:0.01:0.7
    fm = interp1(rho,ff,rad(:,1)-r0);
    c = lsqnonneg([ones(size(rad,1),1) fm],rad(:,2)-min(rad(:,2)));
    err(cnt) = sum((rad(:,2)-min(rad(:,2))-[ones(size(rad,1),1) fm]*c).^2);
    r0v(cnt) = r0;
    cnt = cnt+1;
end
r0 = r0v(err==min(err));
fm = interp1(rho,ff,rad(:,1)-r0);
c = lsqnonneg([ones(size(rad,1),1) fm],rad(:,2)-min(rad(:,2)));
c = [ones(size(rad,1),1) fm]\(rad(:,2)-min(rad(:,2)));
plot(rad(:,1),rad(:,2)-min(rad(:,2)),rad(:,1),[ones(size(rad,1),1) fm]*c)
print -dpng -r300 AzimuthalMode1

rad = azi2; cnt = 1; r0v = []; err = [];
for r0 = 0.4:0.01:0.7
    fm = interp1(rho,ff,rad(:,1)-r0);
    c = lsqnonneg([ones(size(rad,1),1) fm],rad(:,2)-min(rad(:,2)));
    err(cnt) = sum((rad(:,2)-min(rad(:,2))-[ones(size(rad,1),1) fm]*c).^2);
    r0v(cnt) = r0;
    cnt = cnt+1;
end
r0 = r0v(err==min(err));
fm = interp1(rho,ff,rad(:,1)-r0);
c = lsqnonneg([ones(size(rad,1),1) fm],rad(:,2)-min(rad(:,2)));
c = [ones(size(rad,1),1) fm]\(rad(:,2)-min(rad(:,2)));
plot(rad(:,1),rad(:,2)-min(rad(:,2)),rad(:,1),[ones(size(rad,1),1) fm]*c)
print -dpng -r300 AzimuthalMode2