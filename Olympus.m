function [exc, rho, z, fx0, fx2, fz] = Olympus(rhofield, zfield, NA, workdist, n0, n1, lambda, over, focpos)

% function [exc, rho, z, fx0, fx2, fz] = Olympus(rhofield, zfield, NA, workdist, n0, n1, lambda, over, focpos)
% calculates the electric field distribution of a focused laser beam on a glass seurface for a water 
% immersion objective
%
% input arguments:
% rhofield - two-element vector containing the minimum and maximum value of the rho coordinate in mum
% zfield   - two-element vector containing the minimum and maximum value of the z coordinate in mum
% NA       - numerical aperture
% n0       - refractive index of lower medium, usually glass = 1.5
% n1       - refractive index of upper medium, usually water = 1.33
% lambda   - wavelength in mum
% over     - aperture filling factor
% workdist - working distance of the objective, usual something aorund 250 mum
% focpos   - fpcus position in mum
%
% output arguments:
% exc      - excitation intensity distribution
% rho, z   - rho and z coordinate grid
% f**      - field component distributions

if nargin<7 | isempty(lambda)
   lambda = 0.635;
end

k0 = 2*pi/lambda*n0; k1 = 2*pi/lambda*n1;
drho = lambda/25; dz = drho;
rhov = rhofield(1) + (0.5:diff(rhofield)/drho)*drho; 
row = ones(size(rhov));
zv = zfield(1) + (0.5:diff(zfield)/dz)*dz; 
[rho, z] = ndgrid(rhov, zv);
   
chimin = 0;
chimax = asin(NA/n1);

chi = [chimin chimax];
c1 = cos(chi);
s1 = sin(chi);
s0 = n1*s1/n0;
c0 = sqrt(1-s0.^2);
phase = k1*zv'*c1 + (workdist-focpos)*ones(length(zv),1)*(k0*c0-k1*c1);
chiv = 0.5:round((max(max(abs(phase)))-min(min(abs(phase))))/2/pi*25);
length(chiv)
dchi = (chimax-chimin)/(length(chiv)+1);
chiv = chimin + chiv*dchi;

if nargin<8 | isempty(over)
    f2w2 = 0;
else
    f2w2 = -log(1-over)/tan(chimax)^2;
end

if nargin<9
   focpos = 0;
end

h = waitbar(0, 'Focus calculation');
fx0 = zeros(size(rho)); fx2 = fx0;
fz = fx0;
for j=1:length(chiv)
   chi = chiv(j);
   c1 = cos(chi);
   s1 = sin(chi);
   fac = dchi*sin(chi)*sqrt(c1)*exp(-f2w2*tan(chi)^2);
   s0 = n1*s1/n0;
   c0 = sqrt(1-s0^2);
   tp = 2*n0*c0/(n0*c1+n1*c0);
   ts = 2*n0*c0/(n0*c0+n1*c1);
   rp = (n1*c0-n0*c1)/(n1*c0+n0*c1);
   rs = (n0*c0-n1*c1)/(n0*c0+n1*c1);
   phase = k1*c1*zv + (workdist-focpos)*(k0*c0-k1*c1);
   ez = exp(i*phase);
   barg = k1*s1*rhov';
   j0 = fac*besselj(0,barg); j1 = fac*besselj(1,barg); j2 = fac*besselj(2,barg);
   fx0 = fx0 + j0*(tp*c1+ts)*ez; % cos(0*phi)-component
   fx2 = fx2 - j2*(tp*c1-ts)*ez; % cos(2*phi)-component
   fz = fz + 2*i*j1*tp*ez*s1;    % cos(phi)-component
   waitbar(j/length(chiv),h);
end
close(h)

exc = abs(fx0).^2 + abs(fx2).^2 + abs(fz).^2/2;
exc = exc/max(max(exc));

pcolor(rho,z,exc); shading interp; axis image; colorbar

