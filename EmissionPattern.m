function [fxx0, fxx2, fxz, byx0, byx2, byz, rho, z] = EmissionPattern(rho, z, NA, n0, n1, lambda, mag, focpos)

% EmissionPattern calculates the emission pattern of dipoles while imaged with a camera
%
%[fxx0, fxx2, fxz, byx0, byx2, byz, rho, z] = EmissionPattern(rho, z, NA, n0, n1, lambda, mag, focpos)
%
% input arguments:
% r        - two-element vector containing the minimum and maximum value of the rho coordinate
% z        - two-element vector containing the minimum and maximum value of the z coordinate
% NA       - numerical aperture
% n0       - refractive index of lower medium
% n1       - refractive index of upper medium
% lam      - wavelength
% mag      - magnification
% focpos   - focus position along optical axis
%
% output arguments:
% f**      - electric field component distributions
% b**      - magnetic field component distributions
% rho, z   - 2d coordinate fields

ni = 1.0; % it is assumed that imaging medium is air

ki = 2*pi/lambda*ni; k0 = 2*pi/lambda*n0; k1 = 2*pi/lambda*n1;

if length(rho)==2 & length(z)==2
    drho = lambda/25; dz = drho;
    rhov = rho(1) + (0.5:(rho(2)-rho(1))/drho)*drho; 
    zv = z(1) + (0:(z(2)-z(1))/dz)*dz; 
    [rho, z] = ndgrid(rhov, zv);
else
    if length(rho)>1
        drho = rho(2,1)-rho(1,1); 
    else
        drho = lambda/25;
    end
    zv = z(1,:);
end
size(rho)
rhov = mag*(rho(1,1):drho:(rho(end,1)));
zv = mag^2*zv;
row = ones(size(rhov));
    
chimin = 0;
chimax = atan(tan(asin(NA/n0))/mag);

chi = [chimin chimax];
ci = cos(chi);
si = sin(chi);
s0 = sin(atan(si./ci*mag));
s1 = n0*s0/n1;
c0 = max([real(sqrt(1-s0.^2)); 1e-1*ones(1,length(chi))]);
c1 = max([real(sqrt(1-s1.^2)); 1e-1*ones(1,length(chi))]);
col = ones(length(zv),1);

phase = (focpos*(n0/n1*c0./c1)-focpos).*(k0-mag^2*ki*(1-ci)) + focpos*((k1^2-k0^2)/k1./c1);
phase = round(((max(abs(phase))-min(abs(phase))+max(max(ki*ci'*zv)-min(ki*ci'*zv))))/2/pi*25)
if phase<200 phase=200; end
chiv = 0.5:phase;
length(chiv)
dchi = (chimax-chimin)/(length(chiv)+1);
chiv = chimin + chiv*dchi;

if nargin<7
    focpos = 0;
end

h = waitbar(0, 'CEF calculation');
fxx0 = zeros(length(rhov),length(zv)); 
fxx2 = fxx0;
fxy = fxx0;
fxz = fxx0;
fyx = fxx0;
fyy0 = fxx0;
fyy2 = fxx0;
fyz = fxx0;
bxx = fxx0;
bxy0 = fxx0;
bxy2 = fxx0;
byx0 = fxx0;
byx2 = fxx0;
byy = fxx0;
bxz = fxx0;
byz = fxx0;
for j=1:length(chiv)
    chi = chiv(j);
    ci = cos(chi);
    si = sin(chi);
    s0 = sin(atan(si/ci*mag));
    s1 = n0*s0/n1;
    if s1<1
        c0 = sqrt(1-s0^2);
        c1 = sqrt(1-s1^2);
        tp = 2*n1*c1/(n1*c0+n0*c1);
        ts = 2*n1*c1/(n1*c1+n0*c0);
        phase = (focpos*n0/n1*c0/c1-focpos)*(k0-mag^2*ki*(1-ci)) + focpos*(k1^2-k0^2)/k1/c1 + ki*ci*zv;
        ez = exp(i*phase);
        fac = dchi*si*sqrt(ni*ci/n1/c1);
        barg = ki*si*rhov';
        j0 = fac*besselj(0,barg); j1 = fac*besselj(1,barg); j2 = fac*besselj(2,barg);
        
        fxx0 = fxx0 + j0*(tp*c1*ci+ts)*ez; % cos(0*phi)-component
        fxx2 = fxx2 - j2*(tp*c1*ci-ts)*ez; % cos(2*phi)-component
        fxz = fxz - 2*i*j1*tp*ci*s1*ez; % cos(1*phi)-component
        
        byx0 = byx0 + j0*(tp*c1+ts*ci)*ez; % cos(0*phi)-component
        byx2 = byx2 - j2*(tp*c1-ts*ci)*ez; % cos(2*phi)-component
        byz = byz - 2*i*j1*tp*s1*ez;   % cos(1*phi)-component
    end

   waitbar(j/length(chiv));
end
close(h)

%phi=(0:100)*2*pi/100;
%j=1; mm=max(max(real((fxx0(:,j)*ones(size(phi))+fxx2(:,j)*cos(2*phi)).*conj(byx0(:,j)*ones(size(phi))+byx2(:,j)*cos(2*phi))))); figure(gcf); for j=1:100; pcolor(rho(:,1)*cos(phi),rho(:,1)*sin(phi),real((fxx0(:,j)*ones(size(phi))+fxx2(:,j)*cos(2*phi)).*conj(byx0(:,j)*ones(size(phi))+byx2(:,j)*cos(2*phi)))); caxis([0 mm]); axis image; shading interp; drawnow; end
%j=1; mm=max(max(real((fxz(:,j)*cos(phi)).*conj(byz(:,j)*cos(phi))))); figure(gcf); for j=1:100; pcolor(rho(:,1)*cos(phi),rho(:,1)*sin(phi),real((fxz(:,j)*cos(0*phi)).*conj(byz(:,j)*cos(0*phi)))); caxis([0 mm]); axis image; shading interp; drawnow; end