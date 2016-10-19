nn = 25; % half number of pixels along one axis of region of interest on CCD
pixel = 8; % CCD pixel size in mum
NA = 1.4; % NA of objective
focpos = 0.6; % defocusing in mum
n0 = 1.51; % ref. index of objective's immersion medium
n = 1.2; % ref. index of layer with Qdot
n1 = 1; % ref. index of topping layer
d0 = []; 
d = 0.1; % thickness of layer with Qdot in mum
d1 = [];
z = 0; % vertical position of Qdot within its layerin mum
lambda = 0.57; % emission wavelength in mum
mag = 150; % magnification
maxrho = 1.1*nn*sqrt(2)*pixel/mag;

[intx, inty, intz, rho, phi, fxx0, fxx2, fxz, byx0, byx2, byz] = ...
    SEPDipole([0 maxrho], z, NA, n0, n, n1, d0, d, d1, lambda, mag, focpos);

% pcolor(cos(phi)*rho,sin(phi)*rho,intx); axis image; shading interp

be = pi/180*45; % azimuthal angle of dark axis
al = pi/180*45; % polar angle of dark axis
darkaxis = [sin(al)*cos(be), sin(al)*sin(be), cos(al)];
v1 = [-sin(al), cos(al), 0]; 
v2 = cross(darkaxis,v1); v2 = v2/sqrt(sum(v2.^2)); 

om = 0; % third Euler angle of Qdot reference system

tmp = cos(om)*v1 + sin(om)*v2;
v2 = -sin(om)*v1 + cos(om)*v2; % axis of second dipole transition
v1 = tmp; % axis of first dipole transition

[int, ex1, ey1, bx1, by1] = SEPImage(acos(v1*[0 0 1]'),angle(v1(1)+i*v1(2)),nn,pixel,rho,fxx0,fxx2,fxz,byx0,byx2,byz);
[int, ex2, ey2, bx2, by2] = SEPImage(acos(v2*[0 0 1]'),angle(v2(1)+i*v2(2)),nn,pixel,rho,fxx0,fxx2,fxz,byx0,byx2,byz);

kappa = sqrt(0.5); % relative strength of electric/magnetic fields of both dipole transitions
int = real((ex1 + kappa*ex2).*conj(by1 + kappa*by2) - (ey1 + kappa*ey2).*conj(bx1 + kappa*bx2));

mim(int)

