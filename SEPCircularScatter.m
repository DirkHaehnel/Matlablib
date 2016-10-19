function [int, rho, phi, fxx0, fxx2, fxz, byx0, byx2, byz, s0, ezr, ezi, v, pc, ps, psi] = SEPCircularScatter(rho, z, NA, n1, n, n2, d1, d, d2, lambda, mag, focpos)

% SEPCircularScatter calculates the emission pattern of scatterers upon
% circular polarized excitation
%
% [intx, inty, intz, rho, phi, fxx0, fxx2, fxz, byx0, byx2, byzs0, dphase, v, pc, ps, psi] = SEPCircularScatter(rho, z, NA, n1, n, n2, d1, d, d2, lambda, mag, focpos)
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
% s0       - aperture variable
% dphase   - wave front deviation
% v        - angular distribution of radiation for vertical dipole
% pc, ps   - cosine and sine part of angular distribution of radiation for parallel dipole 
% psi      - emission angle for angular distribution of radiation 

maxnum = 1e3;
ni = 1.0; % it is assumed that imaging medium is air

ki = 2*pi/lambda*ni; 

if length(rho)==2
    drho = lambda/50;
    rho = rho(1) + (0:(rho(2)-rho(1))/drho)*drho;
elseif length(rho)>2
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

phase = -n1(1).*c0*focpos;
ez = exp(i*2*pi/lambda*phase);

fac = dchi*si.*sqrt(ci./c0); % aplanatic objective

barg = ki*si'*rho;
j0 = besselj(0,barg); j1 = besselj(1,barg); j2 = besselj(2,barg);

ezi = fac.*ci.*ez;
ezr = fac.*ez;

fxx0 = (ezi.*pc+ezr.*ps)*j0; % cos(0*phi)-component
fxx2 = -(ezi.*pc-ezr.*ps)*j2; % cos(2*phi)-component
%fxz =  -2*i*(ezi.*v)*j1; % cos(1*phi)-component

byx0 = (ezr.*pc+ezi.*ps)*j0; % cos(0*phi)-component
byx2 = -(ezr.*pc-ezi.*ps)*j2; % cos(2*phi)-component
%byz = -2*i*(ezr.*v)*j1; % cos(1*phi)-component

phi = 0:pi/100:2*pi; phi = phi'; col = ones(size(phi));
%int = real((col*fxx0).*conj(col*byx0)+(col*fxx2).*conj(col*byx2)); % 
int = real((col*fxx0+col*fxx2).*conj(col*byx0+col*byx2)); % int for x-dipole along x-pol
%intz = col*real(fxz.*conj(byz)); % int for z-dipole along x-pol



