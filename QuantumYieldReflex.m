function int = QuantumYieldReflex(NA, fd, mag, nv, dv, lamex, w0, focpos);

if NA>0
    maxnum = 3e2;
    ni = 1.0; % it is assumed that imaging medium is air

    ki = 2*pi/lamex*ni;

    rho = [0 3];
    drho = lamex/50;
    rho = rho(1) + (0:(rho(2)-rho(1))/drho)*drho;
    rho = mag*rho; rho = rho(:)';

    chimin = 0;
    chimax = asin(NA/ni/mag);
    dchi = (chimax-chimin)/maxnum;
    chi = chimin + (0.5:maxnum)*dchi;

    ci = cos(chi);
    si = sin(chi);
    s0 = ni/nv(1)*mag*si;
    c0 = sqrt(1-s0.^2);

    % over = [wd*NA, ast, ww]
    rad = NA*fd*si/sin(chimax); % relation between mode angle and entry ray lateral distance
    om = 1./w0^2/2;
    fac = si.*sqrt(ci).*exp(-rad.^2/w0.^2);

    % taking reflection at the lower coverlside interface into account
    [rp rs] = Fresnel(nv(1)*c0,nv,2*pi*dv/lamex);
    barg = ki*si'*rho;

    j0 = besselj(0,barg); j2 = besselj(2,barg);

    ez = exp(-4*pi*i*nv(1)*c0*focpos/lamex);

    fx0 = (fac.*(rp.*ci+rs).*ez)*j0; % cos(0*phi)-component
    fx2 = -(fac.*(rp.*ci-rs).*ez)*j2; % cos(2*phi)-component

    by0 = (fac.*(rs.*ci+rp).*ez)*j0; % cos(0*phi)-component
    by2 = (fac.*(rs.*ci-rp).*ez)*j2; % cos(2*phi)-component

    int = sum(rho.*real(fx0.*conj(by0) + fx2.*conj(by2)));
else
    rp = Fresnel(nv(1),nv,2*pi*dv/lamex);
    int = abs(rp).^2
end