function [intx, intz, rho, phi, fxx0, fxx2, fxz, byx0, byx2, byz] = SEPAniso(rho, z, NA, n1, n, n2, d, lambda, mag, focpos)

% SEPAniso calculates the emission pattern of dipoles within an anisotropic layer while imaged with a camera
%
% [intx, intz, rho, phi, fxx0, fxx2, fxz, byx0, byx2, byz] = SEPAniso(rho, z, NA, n1, n, n2, d1, d, d2, lambda, mag, focpos)
%
% input arguments:
% r        - two-element vector containing the minimum and maximum value of the rho coordinate
% z        - distance of molecule from surface
% NA       - numerical aperture
% n1       - refractive index below molecule
% n        - birefringent refractive index of molecule's layer [no, ne]
% n2       - refractive index above molecule
% d        - thickness of embedding layer
% lambda   - wavelength
% mag      - magnification
% focpos   - focus position along optical axis
%
% output arguments:
% int*     - intensity distribution
% f**      - electric field component distributions
% b**      - magnetic field component distributions
% rho, phi - 2d coordinate fields

maxnum = 1e3;
ni = 1.0; % it is assumed that imaging medium is air

ki = 2*pi/lambda*ni; 

if length(rho)==2
    drho = lambda/50;
    rho = rho(1) + (0:(rho(2)-rho(1))/drho)*drho;
else
    drho = diff(rho(1:2));
end
rho = mag*rho; rho = rho(:)';
row = ones(size(rho));
    
chimin = 0;
chimax = asin(ni*NA/mag);

dchi = (chimax-chimin)/maxnum;
chi = chimin + (0.5:maxnum)*dchi;

ci = cos(chi);
si = sin(chi);
s0 = ni/n1*mag*si;
c0 = sqrt(1-s0.^2);

[v,pc,ps] = DipoleAniso(asin(s0),2*pi*z/lambda,n1,n,n2,2*pi*d/lambda);

v = v.'; pc = pc.'; ps = ps.';

%plot(asin(s0),real(ps),asin(s0),imag(ps),asin(s0),abs(ps));
%v=abs(v); pc=abs(pc); ps = abs(ps);

phase = -n1*c0*focpos;
ez = exp(i*2*pi/lambda*phase);

fac = dchi*si.*sqrt(ni/n1*ci./c0);

barg = ki*si'*rho;
j0 = besselj(0,barg); j1 = besselj(1,barg); j2 = besselj(2,barg);

ezi = fac.*ci.*ez;
ezr = fac.*ez;

fxx0 = (ezi.*pc+ezr.*ps)*j0; % cos(0*phi)-component
fxx2 = -(ezi.*pc-ezr.*ps)*j2; % cos(2*phi)-component
fxz =  -2*i*(ezi.*v)*j1; % cos(1*phi)-component

byx0 = (ezr.*pc+ezi.*ps)*j0; % cos(0*phi)-component
byx2 = -(ezr.*pc-ezi.*ps)*j2; % cos(2*phi)-component
byz = -2*i*(ezr.*v)*j1; % cos(1*phi)-component

phi = 0:pi/100:2*pi; phi = phi'; col = ones(size(phi));
intx = real((col*fxx0+cos(2*phi)*fxx2).*conj(col*byx0+cos(2*phi)*byx2) + (sin(2*phi)*fxx2).*conj(sin(2*phi)*byx2));
intz = col*real(fxz.*conj(byz));

intx = intx/max(max(max(intx)));
intz = intz/max(max(max(intz)));

%for j=1:51 [oil(:,:,j), oilz(:,:,j), rho, phi] = SurfaceEmissionPattern([0 2], 0, 1.4, 1.5, 1.5, 0.67, 66.66, (j-1)/50); end
%for j=1:51 [air(:,:,j), airz(:,:,j), rho, phi] = SurfaceEmissionPattern([0 2], 0, 1.4, 1.5, 1.0, 0.67, 66.66, (j-1)/50); end
%close all; for j=1:51 subplot(121); pcolor(cos(phi)*rho,sin(phi)*rho,air(:,:,j)); colormap gray; shading interp; axis image;  subplot(122); pcolor(cos(phi)*rho,sin(phi)*rho,oil(:,:,j)); colormap gray; shading interp; axis image; pause(0.1); end

%contourf(rho/66.66,(0:50)/50,log10(squeeze(air(1,:,:))'/max(max(max(air)))),-5:0.2:0);