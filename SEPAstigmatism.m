function [intx, inty, intz, rho, phi, fxx, fyx, fyy, fxz, fyz, bxx, byx, bxy, bxz, byz] = SEPAstigmatism(rho, z, NA, n1, n, n2, d1, d, d2, lambda, mag, focpos, atf, gamma)

% SEPAstigmatism calculates the emission pattern of dipoles while imaged with a camera
%
% [intx, inty, intz, rho, phi, fxx, fyx, fyy, fxz, fyz, bxx, byx, bxy, bxz, byz] = SEPAstigmatism(rho, z, NA, n1, n, n2, d1, d, d2, lambda, mag, focpos, atf, gamma)
%
% input arguments:
% rho        - two-element vector containing the minimum and maximum value of the rho coordinate
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
% atf      - cover slide effects
% gamma    - cylindrical lens parameter
%
% output arguments:
% int*     - intensity distributions
% f**      - electric field component distributions
% b**      - magnetic field component distributions
% rho, phi - 2d coordinate fields

maxnum = 1e3;
maxm = 5;
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
2*pi*z/lambda
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
if nargin>13 & ~isempty(gamma)
    gamma = pi*n1(1)/lambda*gamma;    
else
    gamma = 0;
end
om = 0.5*gamma.*s0.^2;
phase = -n1(1)*c0*focpos;
ez = dchi*si.*sqrt(ci./c0).*exp(i*2*pi/lambda*phase-i*om);
ezi = ci.*ez;

barg = ki*si'*rho;
for jj=0:2*maxm eval(['j' int2str(jj) ' = besselj(' int2str(jj) ',barg);']); end
for jj=0:maxm+1 jom(jj+1,:) = besselj(jj,om); end

f0 = 0.5*(pc.*ci+ps); 
f2 = 0.5*(pc.*ci-ps);
g0 = 0.5*(-ps.*ci-pc); 
g2 = 0.5*(-ps.*ci+pc);
fxx(1,:) = (ez.*(f0.*jom(1,:)-i*f2.*jom(2,:)))*j0; % cos(0*phi)
byx(1,:) = -(ez.*(g0.*jom(1,:)+i*g2.*jom(2,:)))*j0; % cos(0*phi)
for j=1:maxm
    eval(['tmp = j' int2str(2*j) ';']);
    fxx(j+1,:) = i^j*(ez.*(i*f2.*(jom(j,:)-jom(j+2,:)) + 2*f0.*jom(j+1,:)))*tmp; % cos(2*j*phi)
    fyx(j,:) = i^(j+1)*(ez.*f2.*(jom(j,:)+jom(j+2,:)))*tmp; % sin(2*j*phi)
    byx(j+1,:) = i^j*(ez.*(i*g2.*(jom(j,:)-jom(j+2,:)) - 2*g0.*jom(j+1,:)))*tmp; % cos(2*j*phi)
    bxx(j,:) = i^(j-1)*(ez.*(g2.*(jom(j,:)+jom(j+2,:))))*tmp; % sin(2*j*phi)
    eval(['tmp = j' int2str(2*j-1) ';']);
    fxz(j,:) = -i^j*(ez.*ci.*v.*(i*jom(j+1,:)-jom(j,:)))*tmp; % cos((2*j-1)*phi)
    byz(j,:) = -i^j*(ez.*v.*(i*jom(j+1,:)-jom(j,:)))*tmp; % cos((2*j-1)*phi)
    fyz(j,:) = i^j*(ez.*ci.*v.*(i*jom(j+1,:)+jom(j,:)))*tmp; % sin((2*j-1)*phi)
    bxz(j,:) = -i^j*(ez.*v.*(i*jom(j+1,:)+jom(j,:)))*tmp; % sin((2*j-1)*phi)
end

phi = 0:pi/100:2*pi; phi = phi'; col = ones(size(phi));
ex = col*fxx(1,:);
by = col*byx(1,:);
ey = 0*ex;
bx = 0*by;
for j=1:maxm
    ex = ex + cos(2*j*phi)*fxx(j+1,:);
    by = by + cos(2*j*phi)*byx(j+1,:);
    ey = ey + sin(2*j*phi)*fyx(j,:);
    bx = bx + sin(2*j*phi)*bxx(j,:);
end
intx = real(ex.*conj(by))-real(ey.*conj(bx)); % int for x-dipole 

fyy(1,:) = (ez.*(f0.*jom(1,:)+i*f2.*jom(2,:)))*j0; % cos(0*phi)
bxy(1,:) = (ez.*(g0.*jom(1,:)-i*g2.*jom(2,:)))*j0; % cos(0*phi)
for j=1:maxm
    eval(['tmp = j' int2str(2*j) ';']);
    fyy(j+1,:) = i^j*(ez.*(-i*f2.*(jom(j,:)-jom(j+2,:)) + 2*f0.*jom(j+1,:)))*tmp; % cos(2*j*phi)
    bxy(j+1,:) = i^j*(ez.*(i*g2.*(jom(j,:)-jom(j+2,:)) + 2*g0.*jom(j+1,:)))*tmp; % cos(2*j*phi)
end

ey = col*fyy(1,:);
bx = col*bxy(1,:);
ex = 0*ey;
by = 0*bx;
for j=1:maxm
    ey = ey + cos(2*j*phi)*fyy(j+1,:);
    bx = bx + cos(2*j*phi)*bxy(j+1,:);
    ex = ex + sin(2*j*phi)*fyx(j,:);
    by = by - sin(2*j*phi)*bxx(j,:);
end
inty = real(ex.*conj(by))-real(ey.*conj(bx)); % int for y-dipole 

ex = 0*ex; ey = ex; bx = ex; by = ex;
for j=1:maxm
    ex = ex + cos((2*j-1)*phi)*fxz(j,:);
    by = by + cos((2*j-1)*phi)*byz(j,:);
    ey = ey + sin((2*j-1)*phi)*fyz(j,:);
    bx = bx + sin((2*j-1)*phi)*bxz(j,:);
end
intz = real(ex.*conj(by)-ey.*conj(bx));

return

% examples 

rhov = [0 1.5]; 
zv = (0:10:300)/1e3;
NA = 1.4; 
n1 = 1.52;
n = 1.33;
n2 = 1.33;
d1 = [];
d = max(zv);
d2 = [];
lambda = 0.57; 
mag = 100;
focpos = 0;
atf = [];
gamma = -0.5;
for j=1:length(zv)
    z = zv(j);
    [intx, inty, intz, rho, phi] = SEPAstigmatism(rhov, z, NA, n1, n, n2, d1, d, d2, lambda, mag, focpos, atf, gamma);
    subplot(221)
    pcolor(cos(phi)*rho,sin(phi)*rho,intx); axis image; shading interp; axis off
    title('horizontal x-dipole')
    subplot(222)
    pcolor(cos(phi)*rho,sin(phi)*rho,inty); axis image; shading interp; axis off
    title('horizontal y-dipole')
    subplot(223)
    pcolor(cos(phi)*rho,sin(phi)*rho,intz); axis image; shading interp; axis off
    title('vertical z-dipole')    
    subplot(224)
    pcolor(cos(phi)*rho,sin(phi)*rho,intx+inty+intz); axis image; shading interp; axis off
    title('isotropic emitter')    
    colormap hot
    suptitle(['\ith\rm = ' mint2str(1e3*z,3,' ') ' nm'])
    drawnow
    eval(['print -dpng tmp' mint2str(j,2)]);
end


focpos = 0.5;
for j=1:length(zv)
    z = zv(j);
    [intx, inty, intz, rho, phi] = SEPDipole(rhov, z, NA, n1, n, n2, d1, d, d2, lambda, mag, focpos);
    subplot(121)
    pcolor(cos(phi)*rho,sin(phi)*rho,intx+inty); axis image; shading interp; axis off
    title('horizontal dipole')
    subplot(122)
    pcolor(cos(phi)*rho,sin(phi)*rho,intz); axis image; shading interp; axis off
    title('vertical dipole')    
    colormap hot
    suptitle(['\ith\rm = ' mint2str(1e3*z,3,' ') ' nm'])
    drawnow
    eval(['print -dpng bla' mint2str(j,2)]);
end

