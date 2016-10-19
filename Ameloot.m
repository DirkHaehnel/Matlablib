% Program for calculating orientational imaging of diamond sample

ng = 1.52;
nw = 1.33;
nd = 2.5;
nsi = 3.9 + i*0.017;
lamex = 0.57;

theta = 0:pi/500:pi/2;
n1 = ng; d1 = [];
n = nw; d = 0.1/lamex*2*pi;
n2 = [nd nsi]; d2 = 1/lamex*2*pi;
z = d;

[v,pc,ps] = DipoleL(theta,z,n1,n,n2,d1,d,d2);
[lvd,lvu,lpd,lpu,qvd,qvu,qpd,qpu,qv,qp] = LifetimeL(z,n1,n,n2,d1,d,d2);

NA = 1.2;
mag = 60;
focpos = 2;
rhov = [0 2];
[intx, inty, intz, rho, phi] = SEPDipole(rhov, z, NA, n1, n, n2, d1, d, d2, lamex, mag, focpos);
subplot(121); pcolor(cos(phi)*rho,sin(phi)*rho,intx+inty); axis image; shading interp;
subplot(122); pcolor(cos(phi)*rho,sin(phi)*rho,intz); axis image; shading interp;
colormap hot