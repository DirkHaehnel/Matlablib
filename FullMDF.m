function [volx, voly, rho, z, fxc, fxs, fyc, fys, fzc, fzs, qv, qp] = FullMDF(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, lamem, mag, av, zpin, atf, kappa, lt, sat, resolution, ring, maxm)

maxnum = 1e3;
i = sqrt(-1);
ni = 1.0; % it is assumed that imaging medium is air

if isempty(d) | d<max(zfield(:))
    d = max(zfield(:));
end

ki = 2*pi/lamem*ni;
k0 = 2*pi/lamem*n0; k = 2*pi/lamem*n; k1 = 2*pi/lamem*n1;

n0all = n0;
nall = n;
n1all = n1;
if length(n0all)==2*(length(d0)+1)
    n0 = n0(1:length(d0)+1);
end
if length(nall)==2
    n = n(1);
end
if length(n1all)==2*(length(d1)+1)
    n1 = n1(1:length(d1)+1);
end

if nargin<22 | isempty(resolution)
    resolution = 20;
end
if nargin<23 | isempty(ring)
    ring = 0;
end
if nargin<24 | isempty(maxm)
    maxm = 2;
end
[fxc, fxs, fyc, fys, fzc, fzs, rho, z] = GaussExc(rhofield, zfield, NA(1), fd, n0, n, n1, d0, d, d1, lamex, over, focpos, atf, resolution, ring(1), maxm);

if length(rho)>1
    drho = rho(2,1)-rho(1,1);
    rhov = mag*(rho(1,1):drho:(rho(end,1)+av/mag));
    uv = mag*rho(:,1);
else
    drho = z(2)-z(1);
    rhov = mag*(rho:drho:(rho+av/mag));
    uv = mag*rho;
end
zv = z(1,:);

if length(ring)>1
    chimin = asin(ring(2)/ni/mag);
else
    chimin = 0;
end
if length(NA)>1
    chimax = asin(NA(2)/ni/mag); % excluding emission above the critical angle
else
    chimax = asin(NA/ni/mag); % excluding emission above the critical angle
end
dchi = chimax/maxnum;

chi = (chimin:dchi:chimax)';
ci = cos(chi);
si = sin(chi);
s0 = ni/n0(1)*mag*si; % Abbe's sine condition
c0 = sqrt(1-s0.^2);

if length(n0all)==2*(length(d0)+1)
    n0 = n0all(length(d0)+2:end);
end
if length(nall)==2
    n = nall(2);
end
if length(n1all)==2*(length(d1)+1)
    n1 = n1all(length(d1)+2:end);
end
d0 = 2*pi/lamem*d0; d = 2*pi/lamem*d; d1 = 2*pi/lamem*d1;
[v,pc,ps] = DipoleL(asin(s0),2*pi*zv/lamem,n0,n,n1,d0,d,d1);
row = ones(size(zv));

% plot((s0*ones(1,size(v,2))).*abs(v).^2,-(c0*ones(1,size(v,2))).*abs(v).^2)

if ~isempty(atf)
    if length(atf)==1 % accounting for reflection losses when using water immersion; atf = ref index of cover slide
        [tmp, tms, tmp, tms] = Fresnel(n0(1)*c0,n0(1),atf);
        [mp, ms, mp, ms] = Fresnel(sqrt(atf^2-n0(1)^2*s0.^2),atf,n0(1));
    else % + aberration for water immersion; atf(2) = thickness mismatch
        [tmp, tms, tmp, tms] = Fresnel(n0(1)*c0,[n0(1) atf(1) atf(1)], 2*pi*atf(2)/lamem);
        [mp, ms, mp, ms] = Fresnel(sqrt(atf(1)^2-n0(1)^2*s0.^2),[atf(1) n0(1) n0(1)], -2*pi*atf(2)/lamem);
    end
    v = ((tmp.*mp).'*row).*v;
    pc = ((tmp.*mp).'*row).*pc;
    ps = ((tms.*ms).'*row).*ps;
end

phase = rem(-k0(1)*c0*focpos(1) + ki*zpin*ci,2*pi);
ez = exp(i*phase);

fac = dchi*si.*sqrt(ci./c0);
barg = ki*rhov'*si';
j0 = besselj(0,barg); j1 = besselj(1,barg); j2 = besselj(2,barg);

ezi = (fac.*ci.*ez)*row;
ezr = (fac.*ez)*row;

f0 = j0*(ezi.*pc+ezr.*ps); % cos(0*phi)-component
f2 = -j2*(ezi.*pc-ezr.*ps); % cos(2*phi)-component
f1 =  -2*i*j1*(ezi.*v); % cos(1*phi)-component

g0 D= j0*(ezr.*pc+ezi.*ps); % cos(0*phi)-component
g2 = -j2*(ezr.*pc-ezi.*ps); % cos(2*phi)-component
g1 = -2*i*j1*(ezr.*v); % cos(1*phi)-component

f00 = real(f0.*conj(g0));
f11 = real(f1.*conj(g1));
f22 = real(f2.*conj(g2));
f01 = real(f0.*conj(g1)+f1.*conj(g0));
f02 = real(f0.*conj(g2)+f2.*conj(g0));
f12 = real(f1.*conj(g2)+f2.*conj(g1));

col = ones(size(uv)); row = ones(size(rhov));
warning off MATLAB:divideByZero
phi = real(acos((av^2-col*rhov.^2-uv.^2*row)./(2*uv*rhov)));
warning on MATLAB:divideByZero
tmp = ones(size(uv))*rhov;
s0 = 2*(pi-phi).*tmp;
s1 = 2*sin(phi).*tmp;
s2 = sin(2*phi).*tmp;
s3 = 2/3*sin(3*phi).*tmp;
s4 = 1/2*sin(4*phi).*tmp;

if isempty(kappa)
    kappa = 0;
end

if isempty(lt) | lt==0
    lam1 = 3/2/n;
    lam2 = 1/2/n;
    lam3 = 3/10/n;
    qv = []; qp = [];
else
    for j=1:length(zv)
        [lvd(j),lvu(j),lpd(j),lpu(j),qvd(j),qvu(j),qpd(j),qpu(j)] = LifetimeL(2*pi*zv(j)/lamem,n0,n,n1,d0,d,d1);
    end
    qv = qvd+qvu;
    qp = qpd+qpu;
    qd = sqrt((qv-qp)./qp);
    aqd = atan(qd);
    if lt==1
        for j=1:length(zv)
            if abs(qd(j))>1e-5 & ~isnan(qd(j))
                lam1(j) = 2*real(aqd(j)/qd(j))/qp(j);
                lam2(j) = 2*real((1-aqd(j)/qd(j))/qd(j)^2)/qp(j);
                lam3(j) = 2*real(((qd(j)^2-3)*qd(j)+3*aqd(j))/qd(j)^5)/3/qp(j);
            elseif isnan(qd(j))
                lam1(j) = 0;
                lam2(j) = 0;
                lam3(j) = 0;
            else
                lam1(j) = 3/2/n;
                lam2(j) = 1/2/n;
                lam3(j) = 3/10/n;
            end
        end
    elseif lt==2 % calculate weighted lifetime
        for j=1:length(zv)
            if abs(qd(j))>1e-5 & ~isnan(qd(j))
                lam1(j) = (1/qv(j)+real(aqd(j)/qd(j))/qp(j))/qp(j);
                lam2(j) = real(aqd(j)/qd(j)^3/qp(j)-1/qd(j)^2/qv(j))/qp(j);
                lam3(j) = real(((qp(j)+2*qv(j))*qd(j) -3*qv(j)*aqd(j))/qd(j)^5)/qp(j)^2/qv(j);
            elseif isnan(qd(j))
                lam1(j) = 0;
                lam2(j) = 0;
                lam3(j) = 0;
            else
                lam1(j) = 9/8/n^2;
                lam2(j) = 3/8/n^2;
                lam3(j) = 9/40/n^2;
            end
        end
    end
    lam1 = diag(lam1); lam2 = diag(lam2); lam3 = diag(lam3);
end

if isempty(sat)
    sat = 0;
end

if sat>0
    tmp = 0; % following procedure not always correct
    for j=0:2*maxm
        psi = j/(2*maxm+1)*2*pi;
        fx = 0*rho; fy = fx; fz = fx;
        for k=0:size(fxc,3)-1
            fx = fx + fxc(:,:,k+1)*cos(k*psi);
            fy = fy + fyc(:,:,k+1)*cos(k*psi);
            fz = fz + fzc(:,:,k+1)*cos(k*psi);
        end
        for k=1:size(fxs,3)
            fx = fx + fxs(:,:,k)*sin(k*psi);
            fy = fy + fys(:,:,k)*sin(k*psi);
            fz = fz + fzs(:,:,k)*sin(k*psi);
        end
        tmp = max([tmp, max(max(abs(fx).^2+abs(fy).^2+abs(fz).^2))]);
    end
    sat = sat/tmp; % normalization of optical saturation to maximum of intensity excitation distribution
end

volx = zeros(size(fxc,1),size(fxc,2),2*maxm+1); voly = volx;
for j=0:2*maxm
    psi = j/(2*maxm+1)*2*pi;
    fx = 0*rho; fy = fx; fz = fx;
    for k=0:maxm
        fx = fx + fxc(:,:,k+1)*cos(k*psi);
        fy = fy + fyc(:,:,k+1)*cos(k*psi);
        fz = fz + fzc(:,:,k+1)*cos(k*psi);
    end
    for k=1:maxm
        fx = fx + fxs(:,:,k)*sin(k*psi);
        fy = fy + fys(:,:,k)*sin(k*psi);
        fz = fz + fzs(:,:,k)*sin(k*psi);
    end

    tmp0 = 1 + sat*(abs(fx).^2+abs(fy).^2+abs(fz).^2);

    u1 = (1/4*(s0*f00).*abs(fx).^2 + 5/16*kappa*(s0*f00).*abs(fx).^2 + 1/4*(s0*f22).*abs(fx).^2 + 1/8*kappa*(s0*f22).*abs(fx).^2 + 1/4*(s0*f00).*abs(fy).^2 - 1/16*kappa*(s0*f00).*abs(fy).^2 + 1/4*(s0*f22).*abs(fy).^2 + 1/8*kappa*(s0*f22).*abs(fy).^2 + 1/4*(s0*f00).*abs(fz).^2 - 1/4*kappa*(s0*f00).*abs(fz).^2 + 1/4*(s0*f22).*abs(fz).^2 - 1/4*kappa*(s0*f22).*abs(fz).^2) + 1/16*(-4*cos(2*psi).*(s2*f02).*abs(fx).^2 - 5*kappa*cos(2*psi).*(s2*f02).*abs(fx).^2 + 3*kappa*cos(4*psi).*(s4*f22).*abs(fx).^2 - 6*kappa*sin(2*psi).*(s2*f02).*real(fx.*conj(fy)) + 6*kappa*sin(4*psi).*(s4*f22).*real(fx.*conj(fy)) - 4*cos(2*psi).*(s2*f02).*abs(fy).^2 + kappa*cos(2*psi).*(s2*f02).*abs(fy).^2 - 3*kappa*cos(4*psi).*(s4*f22).*abs(fy).^2 - 4*cos(2*psi).*(s2*f02).*abs(fz).^2 + 4*kappa*cos(2*psi).*(s2*f02).*abs(fz).^2);
    u2 = 1/4*cos(2*psi).*(s2*f02).*abs(fx).^2 + 7/8*kappa*cos(2*psi).*(s2*f02).*abs(fx).^2 + 1/4*cos(2*psi).*(s2*f11).*abs(fx).^2 + 1/8*kappa*cos(2*psi).*(s2*f11).*abs(fx).^2 - 3/8*kappa*cos(4*psi).*(s4*f22).*abs(fx).^2 + 3/4*kappa*sin(2*psi).*(s2*f02).*real(fx.*conj(fy)) - 3/4*kappa*sin(4*psi).*(s4*f22).*real(fx.*conj(fy)) + 1/4*cos(2*psi).*(s2*f02).*abs(fy).^2 + 1/8*kappa*cos(2*psi).*(s2*f02).*abs(fy).^2 + 1/4*cos(2*psi).*(s2*f11).*abs(fy).^2 + 1/8*kappa*cos(2*psi).*(s2*f11).*abs(fy).^2 + 3/8*kappa*cos(4*psi).*(s4*f22).*abs(fy).^2 + 3/2*kappa*cos(psi).*(s1*f01).*real(fx.*conj(fz)) - 3/4*kappa*cos(psi).*(s1*f12).*real(fx.*conj(fz)) - 3/4*kappa*cos(3*psi).*(s3*f12).*real(fx.*conj(fz)) - 3/4*kappa*sin(psi).*(s1*f12).*real(fy.*conj(fz)) - 3/4*kappa*sin(3*psi).*(s3*f12).*real(fy.*conj(fz)) + 1/4*cos(2*psi).*(s2*f02).*abs(fz).^2 - kappa*cos(2*psi).*(s2*f02).*abs(fz).^2 + 1/4*cos(2*psi).*(s2*f11).*abs(fz).^2 - 1/4*kappa*cos(2*psi).*(s2*f11).*abs(fz).^2 + (-1/4*(s0*f00).*abs(fx).^2 - 7/8*kappa*(s0*f00).*abs(fx).^2 + 1/4*(s0*f11).*abs(fx).^2 + 1/8*kappa*(s0*f11).*abs(fx).^2 - 1/4*(s0*f22).*abs(fx).^2 - 1/2*kappa*(s0*f22).*abs(fx).^2 - 1/4*(s0*f00).*abs(fy).^2 - 1/8*kappa*(s0*f00).*abs(fy).^2 + 1/4*(s0*f11).*abs(fy).^2 + 1/8*kappa*(s0*f11).*abs(fy).^2 - 1/4*(s0*f22).*abs(fy).^2 - 1/2*kappa*(s0*f22).*abs(fy).^2 - 1/4*(s0*f00).*abs(fz).^2 + kappa*(s0*f00).*abs(fz).^2 + 1/4*(s0*f11).*abs(fz).^2 - 1/4*kappa*(s0*f11).*abs(fz).^2 - 1/4*(s0*f22).*abs(fz).^2 + kappa*(s0*f22).*abs(fz).^2);
    u3 = (9/16*kappa*(s0*f00).*abs(fx).^2 - 3/8*kappa*(s0*f11).*abs(fx).^2 + 3/8*kappa*(s0*f22).*abs(fx).^2 + 3/16*kappa*(s0*f00).*abs(fy).^2 - 3/8*kappa*(s0*f11).*abs(fy).^2 + 3/8*kappa*(s0*f22).*abs(fy).^2 - 3/4*kappa*(s0*f00).*abs(fz).^2 + 3/4*kappa*(s0*f11).*abs(fz).^2 - 3/4*kappa*(s0*f22).*abs(fz).^2) - 3/16*(3*kappa*cos(2*psi).*(s2*f02).*abs(fx).^2 + 2*kappa*cos(2*psi).*(s2*f11).*abs(fx).^2 - kappa*cos(4*psi).*(s4*f22).*abs(fx).^2 + 2*kappa*sin(2*psi).*(s2*f02).*real(fx.*conj(fy)) - 2*kappa*sin(4*psi).*(s4*f22).*real(fx.*conj(fy)) + kappa*cos(2*psi).*(s2*f02).*abs(fy).^2 + 2*kappa*cos(2*psi).*(s2*f11).*abs(fy).^2 + kappa*cos(4*psi).*(s4*f22).*abs(fy).^2 + 8*kappa*cos(psi).*(s1*f01).*real(fx.*conj(fz)) - 4*kappa*cos(psi).*(s1*f12).*real(fx.*conj(fz)) - 4*kappa*cos(3*psi).*(s3*f12).*real(fx.*conj(fz)) - 4*kappa*sin(psi).*(s1*f12).*real(fy.*conj(fz)) - 4*kappa*sin(3*psi).*(s3*f12).*real(fy.*conj(fz)) - 4*kappa*cos(2*psi).*(s2*f02).*abs(fz).^2 - 4*kappa*cos(2*psi).*(s2*f11).*abs(fz).^2);

    tmp = (u1*lam1 + u2*lam2 + u3*lam3)./tmp0/(2*maxm+1);
    volx(:,:,1) = volx(:,:,1) + tmp;
    for k = 1:maxm
        volx(:,:,k+1) = volx(:,:,k+1) + 2*tmp*cos(k*psi);
        volx(:,:,maxm+1+k) = volx(:,:,maxm+k+1) + 2*tmp*sin(k*psi);
    end

    u1 = (1/4*(s0*f00).*abs(fx).^2 - 1/16*kappa*(s0*f00).*abs(fx).^2 + 1/4*(s0*f22).*abs(fx).^2 + 1/8*kappa*(s0*f22).*abs(fx).^2 + 1/4*(s0*f00).*abs(fy).^2 + 5/16*kappa*(s0*f00).*abs(fy).^2 + 1/4*(s0*f22).*abs(fy).^2 + 1/8*kappa*(s0*f22).*abs(fy).^2 + 1/4*(s0*f00).*abs(fz).^2 - 1/4*kappa*(s0*f00).*abs(fz).^2 + 1/4*(s0*f22).*abs(fz).^2 - 1/4*kappa*(s0*f22).*abs(fz).^2) + 1/16*(4*cos(2*psi).*(s2*f02).*abs(fx).^2 - kappa*cos(2*psi).*(s2*f02).*abs(fx).^2 - 3*kappa*cos(4*psi).*(s4*f22).*abs(fx).^2 - 6*kappa*sin( 2*psi).*(s2*f02).*real(fx.*conj(fy)) - 6*kappa*sin( 4*psi).*(s4*f22).*real(fx.*conj(fy)) + 4*cos(2*psi).*(s2*f02).*abs(fy).^2 + 5*kappa*cos(2*psi).*(s2*f02).*abs(fy).^2 + 3*kappa*cos(4*psi).*(s4*f22).*abs(fy).^2 + 4*cos(2*psi).*(s2*f02).*abs(fz).^2 - 4*kappa*cos( 2*psi).*(s2*f02).*abs(fz).^2);
    u2 = -1/4*cos(2*psi).*(s2*f02).*abs(fx).^2 - 1/8*kappa*cos(2*psi).*(s2*f02).*abs(fx).^2 - 1/4*cos(2*psi).*(s2*f11).*abs(fx).^2 - 1/8*kappa*cos(2*psi).*(s2*f11).*abs(fx).^2 + 3/8*kappa*cos(4*psi).*(s4*f22).*abs(fx).^2 + 3/4*kappa*sin( 2*psi).*(s2*f02).*real(fx.*conj(fy)) + 3/4*kappa*sin( 4*psi).*(s4*f22).*real(fx.*conj(fy)) - 1/4*cos(2*psi).*(s2*f02).*abs(fy).^2 - 7/8*kappa*cos(2*psi).*(s2*f02).*abs(fy).^2 - 1/4*cos(2*psi).*(s2*f11).*abs(fy).^2 - 1/8*kappa*cos(2*psi).*(s2*f11).*abs(fy).^2 - 3/8*kappa*cos(4*psi).*(s4*f22).*abs(fy).^2 - 3/4*kappa*cos(psi).*(s1*f12).*real(fx.*conj(fz)) + 3/4*kappa*cos( 3*psi).*(s3*f12).*real(fx.*conj(fz)) + 3/2*kappa*sin(psi).*(s1*f01).*real(fy.*conj(fz)) - 3/4*kappa*sin(psi).*(s1*f12).*real(fy.*conj(fz)) + 3/4*kappa*sin( 3*psi).*(s3*f12).*real(fy.*conj(fz)) - 1/4*cos(2*psi).*(s2*f02).*abs(fz).^2 + kappa*cos(2*psi).*(s2*f02).*abs(fz).^2 - 1/4*cos(2*psi).*(s2*f11).*abs(fz).^2 + 1/4*kappa*cos(2*psi).*(s2*f11).*abs(fz).^2 + (-1/4*(s0*f00).*abs(fx).^2 - 1/8*kappa*(s0*f00).*abs(fx).^2 + 1/4*(s0*f11).*abs(fx).^2 + 1/8*kappa*(s0*f11).*abs(fx).^2 - 1/4*(s0*f22).*abs(fx).^2 - 1/2*kappa*(s0*f22).*abs(fx).^2 - 1/4*(s0*f00).*abs(fy).^2 - 7/8*kappa*(s0*f00).*abs(fy).^2 + 1/4*(s0*f11).*abs(fy).^2 + 1/8*kappa*(s0*f11).*abs(fy).^2 - 1/4*(s0*f22).*abs(fy).^2 - 1/2*kappa*(s0*f22).*abs(fy).^2 - 1/4*(s0*f00).*abs(fz).^2 + kappa*(s0*f00).*abs(fz).^2 + 1/4*(s0*f11).*abs(fz).^2 - 1/4*kappa*(s0*f11).*abs(fz).^2 - 1/4*(s0*f22).*abs(fz).^2 + kappa*(s0*f22).*abs(fz).^2);
    u3 = (3/16*kappa*(s0*f00).*abs(fx).^2 - 3/8*kappa*(s0*f11).*abs(fx).^2 + 3/8*kappa*(s0*f22).*abs(fx).^2 + 9/16*kappa*(s0*f00).*abs(fy).^2 - 3/8*kappa*(s0*f11).*abs(fy).^2 + 3/8*kappa*(s0*f22).*abs(fy).^2 - 3/4*kappa*(s0*f00).*abs(fz).^2 + 3/4*kappa*(s0*f11).*abs(fz).^2 - 3/4*kappa*(s0*f22).*abs(fz).^2) + 3/16*(kappa*cos(2*psi).*(s2*f02).*abs(fx).^2 + 2*kappa*cos(2*psi).*(s2*f11).*abs(fx).^2 - kappa*cos(4*psi).*(s4*f22).*abs(fx).^2 - 2*kappa*sin( 2*psi).*(s2*f02).*real(fx.*conj(fy)) - 2*kappa*sin( 4*psi).*(s4*f22).*real(fx.*conj(fy)) + 3*kappa*cos(2*psi).*(s2*f02).*abs(fy).^2 + 2*kappa*cos(2*psi).*(s2*f11).*abs(fy).^2 + kappa*cos(4*psi).*(s4*f22).*abs(fy).^2 + 4*kappa*cos(psi).*(s1*f12).*real(fx.*conj(fz)) - 4*kappa*cos(3*psi).*(s3*f12).*real(fx.*conj(fz)) - 8*kappa*sin(psi).*(s1*f01).*real(fy.*conj(fz)) + 4*kappa*sin(psi).*(s1*f12).*real(fy.*conj(fz)) - 4*kappa*sin( 3*psi).*(s3*f12).*real(fy.*conj(fz)) - 4*kappa*cos(2*psi).*(s2*f02).*abs(fz).^2 - 4*kappa*cos(2*psi).*(s2*f11).*abs(fz).^2);

    tmp = (u1*lam1 + u2*lam2 + u3*lam3)./tmp0/(2*maxm+1);
    voly(:,:,1) = voly(:,:,1) + tmp;
    for k = 1:maxm
        voly(:,:,k+1) = voly(:,:,k+1) + 2*tmp*cos(k*psi);
        voly(:,:,maxm+1+k) = voly(:,:,maxm+k+1) + 2*tmp*sin(k*psi);
    end
end

