load metals

rhofield = [0 0.4];
zfield = 0;
NA = 1.25;
lamex = 0.488;
fd = 164.5e3/100;
over = 4e3;
focpos = 0;
atf = [];
resolution = 100;

% radial polarizartion
n0 = [1.52 silver(wavelength==1e3*lamex)];
n = [1.485 1.54 1.56];
n1 = [silver(wavelength==1e3*lamex) 1];
d0 = [0.04];
d = [0.03 0.04 0.1];
d1 = [0.06];

ring = ['cos(psi)']; maxm = 3;
[fxc, fxs, fyc, fys, fzc, fzs, rho1, z1] = GaussExc(rhofield, [0 d(1)], NA, fd, n0, n(1), [n(2:3) n1], d0, d(1), [d(2:3) d1], lamex, over, focpos, atf, resolution, ring, maxm);
[fxc2, fxs2, fyc2, fys2, fzc2, fzs2] = RotateEMField(fxc, fxs, fyc, fys, fzc, fzs, pi/2);
fx1 = cat(3, fxc + fxc2, fxs + fxs2);
fy1 = cat(3, fyc + fyc2, fys + fys2);
fz1 = cat(3, fzc + fzc2, fzs + fzs2);
[fxc, fxs, fyc, fys, fzc, fzs, rho2, z2] = GaussExc(rhofield, [0 d(2)], NA, fd-d(1), [n0 n(1)], n(2), [n(3) n1], [d0 d(1)], d(2), [d(3) d1], lamex, over, focpos, atf, resolution, ring, maxm);
[fxc2, fxs2, fyc2, fys2, fzc2, fzs2] = RotateEMField(fxc, fxs, fyc, fys, fzc, fzs, pi/2);
fx2 = cat(3, fxc + fxc2, fxs + fxs2);
fy2 = cat(3, fyc + fyc2, fys + fys2);
fz2 = cat(3, fzc + fzc2, fzs + fzs2);
[fxc, fxs, fyc, fys, fzc, fzs, rho3, z3] = GaussExc(rhofield, [0 d(3)], NA, fd-sum(d(1:2)), [n0 n(2:3)], n(3), n1, [d0 d(1:2)], d(3), d1, lamex, over, focpos, atf, resolution, ring, maxm);
[fxc2, fxs2, fyc2, fys2, fzc2, fzs2] = RotateEMField(fxc, fxs, fyc, fys, fzc, fzs, pi/2);
fx3 = cat(3, fxc + fxc2, fxs + fxs2);
fy3 = cat(3, fyc + fyc2, fys + fys2);
fz3 = cat(3, fzc + fzc2, fzs + fzs2);
[fxc, fxs, fyc, fys, fzc, fzs, rho, z] = GaussExc(rhofield, [0 sum(d)], NA, fd, n0, n(3), n1, d0, sum(d), d1, lamex, over, focpos, atf, resolution, ring, maxm);
[fxc2, fxs2, fyc2, fys2, fzc2, fzs2] = RotateEMField(fxc, fxs, fyc, fys, fzc, fzs, pi/2);
fx = cat(3, fxc + fxc2, fxs + fxs2);
fy = cat(3, fyc + fyc2, fys + fys2);
fz = cat(3, fzc + fzc2, fzs + fzs2);


subplot(2,1,1)
FocusImage2D([rho1 rho2 rho3],[z1 d(1)+z2 sum(d(1:2))+z3],[fx1 fx2 fx3])
colorbar
subplot(2,1,2)
FocusImage2D(rho,z,fx)
colorbar
print -dpng -r300 RadialRhoField

subplot(2,1,1)
FocusImage2D([rho1 rho2 rho3],[z1 d(1)+z2 sum(d(1:2))+z3],[fz1 fz2 fz3])
colorbar
subplot(2,1,2)
FocusImage2D(rho,z,fz)
colorbar
print -dpng -r300 RadialZField

ring = ['-sin(psi)']; maxm = 3;
[fxc, fxs, fyc, fys, fzc, fzs, rho1, z1] = GaussExc(rhofield, [0 d(1)], NA, fd, n0, n(1), [n(2:3) n1], d0, d(1), [d(2:3) d1], lamex, over, focpos, atf, resolution, ring, maxm);
[fxc2, fxs2, fyc2, fys2, fzc2, fzs2] = RotateEMField(fxc, fxs, fyc, fys, fzc, fzs, pi/2);
fx1 = cat(3, fxc + fxc2, fxs + fxs2);
fy1 = cat(3, fyc + fyc2, fys + fys2);
fz1 = cat(3, fzc + fzc2, fzs + fzs2);
[fxc, fxs, fyc, fys, fzc, fzs, rho2, z2] = GaussExc(rhofield, [0 d(2)], NA, fd-d(1), [n0 n(1)], n(2), [n(3) n1], [d0 d(1)], d(2), [d(3) d1], lamex, over, focpos, atf, resolution, ring, maxm);
[fxc2, fxs2, fyc2, fys2, fzc2, fzs2] = RotateEMField(fxc, fxs, fyc, fys, fzc, fzs, pi/2);
fx2 = cat(3, fxc + fxc2, fxs + fxs2);
fy2 = cat(3, fyc + fyc2, fys + fys2);
fz2 = cat(3, fzc + fzc2, fzs + fzs2);
[fxc, fxs, fyc, fys, fzc, fzs, rho3, z3] = GaussExc(rhofield, [0 d(3)], NA, fd-sum(d(1:2)), [n0 n(2:3)], n(3), n1, [d0 d(1:2)], d(3), d1, lamex, over, focpos, atf, resolution, ring, maxm);
[fxc2, fxs2, fyc2, fys2, fzc2, fzs2] = RotateEMField(fxc, fxs, fyc, fys, fzc, fzs, pi/2);
fx3 = cat(3, fxc + fxc2, fxs + fxs2);
fy3 = cat(3, fyc + fyc2, fys + fys2);
fz3 = cat(3, fzc + fzc2, fzs + fzs2);
[fxc, fxs, fyc, fys, fzc, fzs, rho, z] = GaussExc(rhofield, [0 sum(d)], NA, fd, n0, n(3), n1, d0, sum(d), d1, lamex, over, focpos, atf, resolution, ring, maxm);
[fxc2, fxs2, fyc2, fys2, fzc2, fzs2] = RotateEMField(fxc, fxs, fyc, fys, fzc, fzs, pi/2);
fx = cat(3, fxc + fxc2, fxs + fxs2);
fy = cat(3, fyc + fyc2, fys + fys2);
fz = cat(3, fzc + fzc2, fzs + fzs2);

subplot(2,1,1)
FocusImage2D([rho1 rho2 rho3],[z1 d(1)+z2 sum(d(1:2))+z3],[fy1 fy2 fy3])
colorbar
subplot(2,1,2)
FocusImage2D(rho,z,fy)
colorbar
print -dpng -r300 AzimuthalPhiField

