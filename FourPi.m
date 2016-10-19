% program estimates the possibility to obtain a tighter localisation of
% light along the optical axis by superposing two different laser foci

NA = 1.2;
n0 = 1.333;
n = n0; n1 =n0;
d0 = []; d =0; d1 = [];
lamex = 0.635;
over = [NA*3e3 0 5e3];
focpos = 10;
atf = [];
resolution = 40*n;
delta = lamex/resolution;
rhofield = [0 2] - delta/2;
zfield =  [8 12] - delta/2;
ring = 0; 

[fx0, fx2, fz, rho, z] = NewExc(rhofield, zfield, NA, n0, n, n1, d0, d, d1, lamex, over, focpos, atf, resolution, ring);
f0 = abs(fx0+fx2).^2+abs(fz).^2; f0 = f0/max(f0(:));
ff = abs(fx0+fx2+conj(fx0+fx2)).^2+abs(fz+conj(fz)).^2; ff = ff/max(ff(:));
mpc(rho,z,ff,f0)

for j=1:40
    shift = j;
    f = abs(fx0(:,1:end-2*shift)+fx2(:,1:end-2*shift)+fx0(:,1+2*shift:end)+fx2(:,1+2*shift:end)).^2+...
        abs(fz(:,1:end-2*shift)+fz(:,1+2*shift:end)).^2;
    f = f/max(f(:));
    mpc(rho(:,1),z(1,1+shift:end-shift,1),ff(:,1+shift:end-shift),f);
    pause; 
end
