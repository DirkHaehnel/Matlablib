function [fx0, fx2, by0, by2, rho, phi, intx, inty] = ConfocalReflex(rho, NA, mag, n0, ng, dg, lamex, over, zim, atf)

maxnum = 3e3;
ni = 1.0; % it is assumed that imaging medium is air

ki = 2*pi/lamex*ni; k0 = 2*pi/lamex*n0; kg = 2*pi/lamex*ng; 

if length(rho)==2
    drho = lamex/50;
    rho = rho(1) + (0:(rho(2)-rho(1))/drho)*drho;
else
    drho = diff(rho(1:2));
end
rho = mag*rho; rho = rho(:)';
row = ones(size(rho));
    
chimin = 0;
chimax = asin(NA/ni/mag);

dchi = (chimax-chimin)/maxnum;
chi = chimin + (0.5:maxnum)*dchi;

ci = cos(chi);
si = sin(chi);
s0 = ni/n0*mag*si;
c0 = sqrt(1-s0.^2);
sg = ni/ng*mag*si;
cg = sqrt(1-sg.^2);

% over = [wd*NA, ast, ww]
rad = over(1)*si/sin(chimax); % relation between mode angle and entry ray lateral distance
zeta = 2*over(2)*pi*over(3)^2/over(1)^2;
w0 = over(3)/sqrt(1+zeta^2);
zeta = [-zeta zeta];
om = (1-i*zeta)./(1+zeta.^2)./w0^2/2;
fac = si.*sqrt(ci).*exp(-rad.^2*sum(om));
om = -i*diff(om)*rad.^2;

% taking reflection at the lower coverlside interface into account
[rp rs tp ts] = Fresnel(ng*cg,ng,n0);
[rpd rsd tpd tsd] = Fresnel(n0*c0,n0,ng);
ez = exp(2*i*kg*dg*cg);
tp = tpd.*rp./(1-rp.^2.*ez).*tp + rpd./ez;
ts = tsd.*rs./(1-rs.^2.*ez).*ts + rsd./ez;
ez = exp(-2*i*k0*zim*c0);

if nargin>9 & ~isempty(atf)
    ez = ez.*exp(2*i*(kg*cg-k0*c0)*atf);
end
barg = ki*si'*rho;

ez = ez.*exp(i*5*pi*(c0.^2-c0.^3));

if over(2)==0
    j0 = besselj(0,barg); j2 = besselj(2,barg);
  
    fx0 = (fac.*(tp.*ci+ts).*ez)*j0; % cos(0*phi)-component
    fx2 = -(fac.*(tp.*ci-ts).*ez)*j2; % cos(2*phi)-component

    by0 = (fac.*(ts.*ci+tp).*ez)*j0; % cos(0*phi)-component
    by2 = (fac.*(ts.*ci-tp).*ez)*j2; % cos(2*phi)-component
    
    if nargout>5
        phi = 0:pi/100:2*pi; phi = phi'; col = ones(size(phi));
        intx = real((col*fx0+cos(2*phi)*fx2).*conj(col*by0+cos(2*phi)*by2));
        inty = real((sin(2*phi)*fx2).*conj(sin(2*phi)*by2));
    end
else
    for jj=0:2:8 eval(['j' int2str(jj) ' = besselj(' int2str(jj) ',barg);']); end
    for jj=0:4 eval(['jom' int2str(jj) ' = besselj(' int2str(jj) ',om);']); end
    f0 = fac.*(tp.*ci+ts).*ez; 
    f2 = fac.*(tp.*ci-ts).*ez;

    fx0(1,:) = (f0.*jom0 - i*f2.*jom1)*j0; % cos(0*phi)
    fx0(2,:) = (-f2.*jom0 + 2*i*f0.*jom1 + f2.*jom2)*j2; % cos(2*phi)
    fx0(3,:) = (-i*f2.*jom1 - 2*f0.*jom2 + i*f2.*jom3)*j4; % cos(4*phi)
    fx0(4,:) = (f2.*jom2 - 2*i*f0.*jom3 - f2.*jom4)*j6; % cos(6*phi)
    
    fx2(1,:) = -(f2.*(jom0+jom2))*j2; % sin(2*phi)
    fx2(2,:) = -i*(f2.*(jom1+jom3))*j4; % sin(4*phi)
    fx2(3,:) = (f2.*(jom2+jom4))*j6; % sin(6*phi)
    fx2(4,:) = i*(f2.*jom3)*j8; % sin(8*phi)

    f0 = fac.*(ts.*ci+tp).*ez; 
    f2 = fac.*(ts.*ci-tp).*ez;

    by0(1,:) = (f0.*jom0 + i*f2.*jom1)*j0; % cos(0*phi)
    by0(2,:) = (-f2.*jom2 + 2*i*f0.*jom1 + f2.*jom0)*j2; % cos(2*phi)
    by0(3,:) = (-i*f2.*jom3 - 2*f0.*jom2 + i*f2.*jom1)*j4; % cos(4*phi)
    by0(4,:) = (f2.*jom4 - 2*i*f0.*jom3 - f2.*jom2)*j6; % cos(6*phi)
    
    by2(1,:) = -(f2.*(jom0+jom2))*j2; % sin(2*phi)
    by2(2,:) = -i*(f2.*(jom1+jom3))*j4; % sin(4*phi)
    by2(3,:) = (f2.*(jom2+jom4))*j6; % sin(6*phi)
    by2(4,:) = i*(f2.*jom3)*j8; % sin(8*phi)

    if nargout>5
        phi = 0:pi/100:2*pi; phi = phi'; col = ones(size(phi));
        ex = col*fx0(1,:);
        ey = sin(2*phi)*fx2(1,:);
        by = col*by0(1,:);
        bx = sin(2*phi)*by2(1,:);
        for j=1:3
            ex = ex + cos(2*j*phi)*fx0(j+1,:);
            ey = ey + sin(2*(j+1)*phi)*fx2(j+1,:);
            by = by + cos(2*j*phi)*by0(j+1,:);
            bx = bx + sin(2*(j+1)*phi)*by2(j+1,:);
        end
        intx = real(ex.*conj(by));
        inty = -real(ey.*conj(bx));
    end
end

% pcolor(cos(phi)*rho,sin(phi)*rho,intx+inty); shading interp; axis image; colormap gray; axis off
