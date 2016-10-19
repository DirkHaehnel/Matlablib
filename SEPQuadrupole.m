function [intx, inty, intz, rho, phi, fxx0, fxx2, fxz, byx0, byx2, byz, v, pc, ps] = SEPQuadrupole(rho, z, NA, n1, n, n2, d1, d, d2, lambda, mag, focpos, atf, ord)

% SEPQuadrupole calculates the emission pattern of dipoles while imaged with a camera
%
% [intx, intz, rho, phi, fxx0, fxx2, fxz, byx0, byx2, byz] = SEPQuadrupole(rho, z, NA, n1, n2, lambda, mag, focpos)
%
% input arguments:
% r        - two-element vector containing the minimum and maximum value of the rho coordinate
% z        - distance of molecule from surface
% NA       - numerical aperture
% n1       - refractive index vector (bottom to top) below molecule
% n2       - refractive index vector (bottom to top) above molecule
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
chimax = asin(NA/mag/ni);
dchi = (chimax-chimin)/maxnum;
chi = chimin + (0.5:maxnum)*dchi;

ci = cos(chi);
si = sin(chi);
s0 = ni/n1(1)*mag*si;
c0 = sqrt(1-s0.^2);
psi = asin(s0);
if length(ord)==1
    [v,pc,ps] = MultipoleL(psi,2*pi*z/lambda,n1,n,n2,2*pi*d1/lambda,2*pi*d/lambda,2*pi*d2/lambda,ord);
else
    [v,pc,ps] = MultipoleL(psi,2*pi*z/lambda,n1,n,n2,2*pi*d1/lambda,2*pi*d/lambda,2*pi*d2/lambda,0);
    v = ord(1)*v; pc = ord(1)*pc; ps = ord(1)*ps;
    for j=2:length(ord)
        [v1,pc1,ps1] = MultipoleL(psi,2*pi*z/lambda,n1,n,n2,2*pi*d1/lambda,2*pi*d/lambda,2*pi*d2/lambda,j-1);
        v = v + ord(j)*v1;
        pc = pc + ord(j)*pc1;
        ps = ps + ord(j)*ps1;
    end
end
v = v.'; pc = pc.'; ps = ps.';

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

phase = -n1(1)*c0*focpos;
ez = exp(i*2*pi/lambda*phase);

fac = si.*sqrt(ci./c0);

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