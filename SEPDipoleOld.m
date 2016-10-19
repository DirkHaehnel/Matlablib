function [intx, inty, intz, rho, phi, fxx0, fxx2, fxz, byx0, byx2, byz, v, pc, ps, psi] = SEPDipole(rho, z, NA, n1, n, n2, d1, d, d2, lambda, mag, focpos, atf)

% EmissionPattern calculates the emission pattern of dipoles while imaged with a camera
%
% [intx, inty, intz, rho, phi, fxx0, fxx2, fxz, byx0, byx2, byz] = SEPDipole(rho, z, NA, n1, n, n2, d1, d, d2, lambda, mag, focpos)
%
% input arguments:
% r        - two-element vector containing the minimum and maximum value of the rho coordinate
% z        - distance of molecule from surface
% NA       - numerical aperture
% n1       - refractive index vector (bottom to top) below molecule
% n        - refractive index of molecule's layer
% n2       - refractive index vector (bottom to top) above molecule
% d1       - thickness vector for layers below molecule
% d        - thickness of embedding layer
% d2       - thickness vector for layers above molecule
% lambda   - wavelength
% mag      - magnification
% focpos   - focus position along optical axis
%
% output arguments:
% int*     - intensity distributions
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
chimax = asin(NA/mag/ni);
dchi = (chimax-chimin)/maxnum;
chi = chimin + (0.5:maxnum)*dchi;

ci = cos(chi);
si = sin(chi);
s0 = ni/n1(1)*mag*si;
c0 = sqrt(1-s0.^2);
psi = asin(s0);
[v,pc,ps] = DipoleL(psi,2*pi*z/lambda,n1,n,n2,2*pi*d1/lambda,2*pi*d/lambda,2*pi*d2/lambda);
v = v.'; ps = ps.'; pc = pc.';

if nargin>12 & ~isempty(atf)
    if length(atf)==1 % accounting for reflection losses when using water immersion; atf = ref index of cover slide
        [tmp, tms, tmp, tms] = fresnel(n1(1)*c0,n1(1),atf);
        [mp, ms, mp, ms] = fresnel(sqrt(atf^2-n1(1)^2*s0.^2),atf,n1(1));
    else % + aberration for water immersion; atf(2) = thickness mismatch
        [tmp, tms, tmp, tms] = fresnel(n1(1)*c0,[n1(1) atf(1) atf(1)], 2*pi*atf(2)/lambda);        
        [mp, ms, mp, ms] = fresnel(sqrt(atf(1)^2-n1(1)^2*s0.^2),[atf(1) n1(1) n1(1)], -2*pi*atf(2)/lambda);        
    end
    v = tmp.*mp.*v;
    pc = tmp.*mp.*pc;    
    ps = tms.*ms.*ps;
end

phase = -n1(1).*c0*focpos; %-kappa(1)*s0.^2-kappa(2)*s0.^4;
%alternative formulation but same result:
%phase = ni*mag*focpos*s0./c0.*si-n1(1)./c0*focpos; 
ez = exp(i*2*pi/lambda*phase);

fac = si.*sqrt(ci./c0); % aplanatic objective

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
intx = real((col*fxx0+cos(2*phi)*fxx2).*conj(col*byx0+cos(2*phi)*byx2)); % int for x-dipole along x-pol
inty = real((sin(2*phi)*fxx2).*conj(sin(2*phi)*byx2)); % int for x-dipole along y-pol
intz = col*real(fxz.*conj(byz)); % int for z-dipole along x-pol

%int = real((cos(al)*(col*fxx0+cos(2*phi)*fxx2)+sin(al)*cos(phi)*fxz).*conj(cos(al)*(col*byx0+cos(2*phi)*byx2)+sin(al)*cos(phi)*byz) + (cos(al)*sin(2*phi)*fxx2+sin(al)*sin(phi)*fxz).*conj(cos(al)*sin(2*phi)*byx2+sin(al)*sin(phi)*byz));

%for j=1:51 [oilx(:,:,j), oily(:,:,j), oilz(:,:,j), rho, phi] = SEPDipole((0:20)*6.5/60, 0, 1.4, 1.5, 1.5, 1.0, [], 0, [], 0.67, 60, (j-1)/50); end
%for j=1:51 [airx(:,:,j), airy(:,:,j), airz(:,:,j), rho, phi] = SEPDipole([0 2], 0, 1.4, 1.5, 1.0, 1.0, [], 0, [], 0.67, 60, (j-1)/50); end
%close all; for j=1:51 subplot(121); pcolor(cos(phi)*rho,sin(phi)*rho,airx(:,:,j)+airy(:,:,j)); colormap gray; shading interp; axis image;  subplot(122); pcolor(cos(phi)*rho,sin(phi)*rho,oilx(:,:,j)+oily(:,:,j)); colormap gray; shading interp; axis image; pause(0.1); end

%close all; for j=1:51 pcolor(cos(phi)*rho,sin(phi)*rho,oilx(:,:,j)+oily(:,:,j)); colormap gray; shading interp; axis image; title(num2str(j)); pause(0.1); end

%close all; for j=2:51 subplot(10,5,j-1); pcolor(cos(phi)*rho,sin(phi)*rho,oilx(:,:,j)+oily(:,:,j)); axis off; colormap gray; shading interp; axis image; end

%contourf(rho/66.66,(0:50)/50,log10(squeeze(airx(1,:,:)+airy(1,:,:))'/max(max(max(airx+airy)))),-5:0.2:0);