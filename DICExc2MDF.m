function mdf = DICExc2MDF(exc, NA, n0, n, n1, focpos, lamem, mag, av, zpin, atf, kappa, lt, pulse, sat, triplet)
 
% important: DIC prism is always assumed to work along x-direction !!! 
% i.e. focus positions have to form a line along x!!!
% first focus is x-polarized, second focus is y-polarized

maxnum = 1e3;
ni = 1.0; % it is assumed that imaging medium is air

ki = 2*pi/lamem*ni;
k0 = 2*pi/lamem*n0; 

if nargin<10
    zpin = 0;
end

if nargin<11 
    atf = [];
end

if nargin<12 || isempty(kappa)
    kappa = 0;
end

if nargin<13 || isempty(lt)
    lt = 0;
end

if nargin<14 || isempty(pulse)
    pulse = [1 1];
end
if length(pulse)==1
    pulse = [pulse pulse];
end

if nargin<15 || isempty(sat)
    sat = 0;
end

if nargin<16 || isempty(triplet)
    triplet = 0;
end

maxm = exc.maxm;
rho = exc.rho;
z = exc.z;
fxc1 = exc.fxc1;
fxs1 = exc.fxs1;
fyc1 = exc.fyc1;
fys1 = exc.fys1;
fzc1 = exc.fzc1;
fzs1 = exc.fzs1;
fxc2 = exc.fxc2;
fxs2 = exc.fxs2;
fyc2 = exc.fyc2;
fys2 = exc.fys2;
fzc2 = exc.fzc2;
fzs2 = exc.fzs2;
b = mag*(exc.focpos(1,2)-exc.focpos(2,2))/2;

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

chimin = 0;
chimax = asin(NA/ni/mag); % excluding emission above the critical angle
dchi = chimax/maxnum;
chi = (chimin:dchi:chimax)';
ci = cos(chi);
si = sin(chi);
s0 = ni/n0(1)*mag*si; % Abbe's sine condition
c0 = sqrt(1-s0.^2);

d0 = 2*pi/lamem*exc.d0; d = 2*pi/lamem*exc.d; d1 = 2*pi/lamem*exc.d1;
[v,pc,ps] = DipoleL(asin(s0),2*pi*zv/lamem,n0,n,n1,d0,d,d1);
row = ones(size(zv));

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
ez = exp(1i*phase);

fac = dchi*si.*sqrt(ci./c0);
barg = ki*rhov'*si';
j0 = besselj(0,barg); j1 = besselj(1,barg); j2 = besselj(2,barg);

ezi = (fac.*ci.*ez)*row;
ezr = (fac.*ez)*row;

f0 = j0*(ezi.*pc+ezr.*ps); % cos(0*phi)-component
f2 = -j2*(ezi.*pc-ezr.*ps); % cos(2*phi)-component
f1 =  -2*1i*j1*(ezi.*v); % cos(1*phi)-component

g0 = j0*(ezr.*pc+ezi.*ps); % cos(0*phi)-component
g2 = -j2*(ezr.*pc-ezi.*ps); % cos(2*phi)-component
g1 = -2*1i*j1*(ezr.*v); % cos(1*phi)-component

f00 = real(f0.*conj(g0));
f11 = real(f1.*conj(g1));
f22 = real(f2.*conj(g2));
f01 = real(f0.*conj(g1)+f1.*conj(g0));
f02 = real(f0.*conj(g2)+f2.*conj(g0));
f12 = real(f1.*conj(g2)+f2.*conj(g1));

warning off MATLAB:divideByZero
warning on MATLAB:divideByZero

if lt==0
    lam1 = 3/2/n;
    lam2 = 1/2/n;
    lam3 = 3/10/n;
    qv = []; qp = [];
else
    for j=1:length(zv)
        [tmp,tmp,tmp,tmp,qvd(j),qvu(j),qpd(j),qpu(j)] = LifetimeL(2*pi*zv(j)/lamem,n0,n,n1,d0,d,d1);
    end
    qv = qvd+qvu;
    qp = qpd+qpu;
    qd = sqrt((qv-qp)./qp);
    aqd = atan(qd);
    if lt==1
        for j=1:length(zv)
            if abs(qd(j))>1e-5 && ~isnan(qd(j))
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
            if abs(qd(j))>1e-5 && ~isnan(qd(j))
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

if sat>0
    tmp = 0; % following procedure is only approximative
    for j=0:maxm
        psi = j/(2*maxm+1)*2*pi;
        fx = 0*rho; fy = fx; fz = fx;
        for k=0:size(fxc1,3)-1
            fx = fx + fxc1(:,:,k+1)*cos(k*psi);
            fy = fy + fyc1(:,:,k+1)*cos(k*psi);
            fz = fz + fzc1(:,:,k+1)*cos(k*psi);
        end
        for k=1:size(fxs1,3)
            fx = fx + fxs1(:,:,k)*sin(k*psi);
            fy = fy + fys1(:,:,k)*sin(k*psi);
            fz = fz + fzs1(:,:,k)*sin(k*psi);
        end
        tmp = max([tmp, max(max(abs(fx).^2+abs(fy).^2+abs(fz).^2))]);
    end
    sat1 = sat/tmp;
    tmp = 0;
    for j=0:maxm
        psi = j/(2*maxm+1)*2*pi;
        fx = 0*rho; fy = fx; fz = fx;
        for k=0:size(fxc2,3)-1
            fx = fx + fxc2(:,:,k+1)*cos(k*psi);
            fy = fy + fyc2(:,:,k+1)*cos(k*psi);
            fz = fz + fzc2(:,:,k+1)*cos(k*psi);
        end
        for k=1:size(fxs2,3)
            fx = fx + fxs2(:,:,k)*sin(k*psi);
            fy = fy + fys2(:,:,k)*sin(k*psi);
            fz = fz + fzs2(:,:,k)*sin(k*psi);
        end
        tmp = max([tmp, max(max(abs(fx).^2+abs(fy).^2+abs(fz).^2))]);
    end
    sat2 = sat/tmp;
else
    sat1 = 0; sat2 = 0;
end

col = ones(size(uv)); 
row = ones(size(rhov)); 
mat = ones(size(uv))*rhov; 
warning off MATLAB:divideByZero
volx1 = zeros(size(rho,1),size(rho,2),2*maxm+1); 
voly1 = volx1; volx2 = volx1; voly2 = volx1;
zrow = ones(1,size(rho,2));
for j=0:4*maxm
    psi = j/(4*maxm+1)*2*pi;

    % first focus
    fx = fxc1(:,:,1); fy = fyc1(:,:,1); fz = fzc1(:,:,1);
    for k=1:maxm
        fx = fx + fxc1(:,:,k+1)*cos(k*psi) + fxs1(:,:,k)*sin(k*psi);
        fy = fy + fyc1(:,:,k+1)*cos(k*psi) + fys1(:,:,k)*sin(k*psi);
        fz = fz + fzc1(:,:,k+1)*cos(k*psi) + fzs1(:,:,k)*sin(k*psi);
    end

    if sat1>0
        int = sat1*(abs(fx).^2+abs(fy).^2+abs(fz).^2);
        tmp = PulsedExcitation(int,1,pulse(1),pulse(2));
        tmp0 = 1 + triplet*tmp;
        tmp = sqrt(tmp./int);
        fx = fx.*tmp;
        fy = fy.*tmp;
        fz = fz.*tmp;
    else
        tmp0 = 1 + triplet*(abs(fx).^2+abs(fy).^2+abs(fz).^2);
    end

    tmp = sqrt(uv.^2+b.^2-2*b*uv*cos(psi));
    phi = real(acos((av^2-col*rhov.^2-tmp.^2*row)./(2*tmp*rhov)));
    s0 = 2*(pi-phi).*mat;
    s1 = 2*sin(phi).*mat;
    s2 = sin(2*phi).*mat;
    s3 = 2/3*sin(3*phi).*mat;
    s4 = 1/2*sin(4*phi).*mat;

    psidpsi = psi + real(acos((uv-b*cos(psi))./sqrt(uv.^2+b^2-2*b*uv*cos(psi))));
    psidpsi(~isfinite(psidpsi)) = 0;
    ux1 = (1/4*(s0*f00).*abs(fx).^2 + 5/16*kappa*(s0*f00).*abs(fx).^2 + 1/4*(s0*f22).*abs(fx).^2 + 1/8*kappa*(s0*f22).*abs(fx).^2 + 1/4*(s0*f00).*abs(fy).^2 - 1/16*kappa*(s0*f00).*abs(fy).^2 + 1/4*(s0*f22).*abs(fy).^2 + 1/8*kappa*(s0*f22).*abs(fy).^2 + 1/4*(s0*f00).*abs(fz).^2 - 1/4*kappa*(s0*f00).*abs(fz).^2 + 1/4*(s0*f22).*abs(fz).^2 - 1/4*kappa*(s0*f22).*abs(fz).^2) + 1/16*(-4*(cos(2*psidpsi)*zrow).*(s2*f02).*abs(fx).^2 - 5*kappa*(cos(2*psidpsi)*zrow).*(s2*f02).*abs(fx).^2 + 3*kappa*(cos(4*psidpsi)*zrow).*(s4*f22).*abs(fx).^2 - 6*kappa*(sin(2*psidpsi)*zrow).*(s2*f02).*real(fx.*conj(fy)) + 6*kappa*(sin(4*psidpsi)*zrow).*(s4*f22).*real(fx.*conj(fy)) - 4*(cos(2*psidpsi)*zrow).*(s2*f02).*abs(fy).^2 + kappa*(cos(2*psidpsi)*zrow).*(s2*f02).*abs(fy).^2 - 3*kappa*(cos(4*psidpsi)*zrow).*(s4*f22).*abs(fy).^2 - 4*(cos(2*psidpsi)*zrow).*(s2*f02).*abs(fz).^2 + 4*kappa*(cos(2*psidpsi)*zrow).*(s2*f02).*abs(fz).^2);
    ux2 = 1/4*(cos(2*psidpsi)*zrow).*(s2*f02).*abs(fx).^2 + 7/8*kappa*(cos(2*psidpsi)*zrow).*(s2*f02).*abs(fx).^2 + 1/4*(cos(2*psidpsi)*zrow).*(s2*f11).*abs(fx).^2 + 1/8*kappa*(cos(2*psidpsi)*zrow).*(s2*f11).*abs(fx).^2 - 3/8*kappa*(cos(4*psidpsi)*zrow).*(s4*f22).*abs(fx).^2 + 3/4*kappa*(sin(2*psidpsi)*zrow).*(s2*f02).*real(fx.*conj(fy)) - 3/4*kappa*(sin(4*psidpsi)*zrow).*(s4*f22).*real(fx.*conj(fy)) + 1/4*(cos(2*psidpsi)*zrow).*(s2*f02).*abs(fy).^2 + 1/8*kappa*(cos(2*psidpsi)*zrow).*(s2*f02).*abs(fy).^2 + 1/4*(cos(2*psidpsi)*zrow).*(s2*f11).*abs(fy).^2 + 1/8*kappa*(cos(2*psidpsi)*zrow).*(s2*f11).*abs(fy).^2 + 3/8*kappa*(cos(4*psidpsi)*zrow).*(s4*f22).*abs(fy).^2 + 3/2*kappa*(cos(psidpsi)*zrow).*(s1*f01).*real(fx.*conj(fz)) - 3/4*kappa*(cos(psidpsi)*zrow).*(s1*f12).*real(fx.*conj(fz)) - 3/4*kappa*(cos(3*psidpsi)*zrow).*(s3*f12).*real(fx.*conj(fz)) - 3/4*kappa*(sin(psidpsi)*zrow).*(s1*f12).*real(fy.*conj(fz)) - 3/4*kappa*(sin(3*psidpsi)*zrow).*(s3*f12).*real(fy.*conj(fz)) + 1/4*(cos(2*psidpsi)*zrow).*(s2*f02).*abs(fz).^2 - kappa*(cos(2*psidpsi)*zrow).*(s2*f02).*abs(fz).^2 + 1/4*(cos(2*psidpsi)*zrow).*(s2*f11).*abs(fz).^2 - 1/4*kappa*(cos(2*psidpsi)*zrow).*(s2*f11).*abs(fz).^2 + (-1/4*(s0*f00).*abs(fx).^2 - 7/8*kappa*(s0*f00).*abs(fx).^2 + 1/4*(s0*f11).*abs(fx).^2 + 1/8*kappa*(s0*f11).*abs(fx).^2 - 1/4*(s0*f22).*abs(fx).^2 - 1/2*kappa*(s0*f22).*abs(fx).^2 - 1/4*(s0*f00).*abs(fy).^2 - 1/8*kappa*(s0*f00).*abs(fy).^2 + 1/4*(s0*f11).*abs(fy).^2 + 1/8*kappa*(s0*f11).*abs(fy).^2 - 1/4*(s0*f22).*abs(fy).^2 - 1/2*kappa*(s0*f22).*abs(fy).^2 - 1/4*(s0*f00).*abs(fz).^2 + kappa*(s0*f00).*abs(fz).^2 + 1/4*(s0*f11).*abs(fz).^2 - 1/4*kappa*(s0*f11).*abs(fz).^2 - 1/4*(s0*f22).*abs(fz).^2 + kappa*(s0*f22).*abs(fz).^2);
    ux3 = (9/16*kappa*(s0*f00).*abs(fx).^2 - 3/8*kappa*(s0*f11).*abs(fx).^2 + 3/8*kappa*(s0*f22).*abs(fx).^2 + 3/16*kappa*(s0*f00).*abs(fy).^2 - 3/8*kappa*(s0*f11).*abs(fy).^2 + 3/8*kappa*(s0*f22).*abs(fy).^2 - 3/4*kappa*(s0*f00).*abs(fz).^2 + 3/4*kappa*(s0*f11).*abs(fz).^2 - 3/4*kappa*(s0*f22).*abs(fz).^2) - 3/16*(3*kappa*(cos(2*psidpsi)*zrow).*(s2*f02).*abs(fx).^2 + 2*kappa*(cos(2*psidpsi)*zrow).*(s2*f11).*abs(fx).^2 - kappa*(cos(4*psidpsi)*zrow).*(s4*f22).*abs(fx).^2 + 2*kappa*(sin(2*psidpsi)*zrow).*(s2*f02).*real(fx.*conj(fy)) - 2*kappa*(sin(4*psidpsi)*zrow).*(s4*f22).*real(fx.*conj(fy)) + kappa*(cos(2*psidpsi)*zrow).*(s2*f02).*abs(fy).^2 + 2*kappa*(cos(2*psidpsi)*zrow).*(s2*f11).*abs(fy).^2 + kappa*(cos(4*psidpsi)*zrow).*(s4*f22).*abs(fy).^2 + 8*kappa*(cos(psidpsi)*zrow).*(s1*f01).*real(fx.*conj(fz)) - 4*kappa*(cos(psidpsi)*zrow).*(s1*f12).*real(fx.*conj(fz)) - 4*kappa*(cos(3*psidpsi)*zrow).*(s3*f12).*real(fx.*conj(fz)) - 4*kappa*(sin(psidpsi)*zrow).*(s1*f12).*real(fy.*conj(fz)) - 4*kappa*(sin(3*psidpsi)*zrow).*(s3*f12).*real(fy.*conj(fz)) - 4*kappa*(cos(2*psidpsi)*zrow).*(s2*f02).*abs(fz).^2 - 4*kappa*(cos(2*psidpsi)*zrow).*(s2*f11).*abs(fz).^2);

    tmp = sqrt(uv.^2+4*b.^2-4*b*uv*cos(psi));
    phi = real(acos((av^2-col*rhov.^2-tmp.^2*row)./(2*tmp*rhov)));
    s0 = 2*(pi-phi).*mat;
    s1 = 2*sin(phi).*mat;
    s2 = sin(2*phi).*mat;
    s3 = 2/3*sin(3*phi).*mat;
    s4 = 1/2*sin(4*phi).*mat;

    psidpsi = psi + real(acos((uv-2*b*cos(psi))./sqrt(uv.^2+4*b^2-4*b*uv*cos(psi))));
    psidpsi(~isfinite(psidpsi)) = 0;
    uy1 = (1/4*(s0*f00).*abs(fx).^2 - 1/16*kappa*(s0*f00).*abs(fx).^2 + 1/4*(s0*f22).*abs(fx).^2 + 1/8*kappa*(s0*f22).*abs(fx).^2 + 1/4*(s0*f00).*abs(fy).^2 + 5/16*kappa*(s0*f00).*abs(fy).^2 + 1/4*(s0*f22).*abs(fy).^2 + 1/8*kappa*(s0*f22).*abs(fy).^2 + 1/4*(s0*f00).*abs(fz).^2 - 1/4*kappa*(s0*f00).*abs(fz).^2 + 1/4*(s0*f22).*abs(fz).^2 - 1/4*kappa*(s0*f22).*abs(fz).^2) + 1/16*(4*(cos(2*psidpsi)*zrow).*(s2*f02).*abs(fx).^2 - kappa*(cos(2*psidpsi)*zrow).*(s2*f02).*abs(fx).^2 - 3*kappa*(cos(4*psidpsi)*zrow).*(s4*f22).*abs(fx).^2 - 6*kappa*(sin( 2*psidpsi)*zrow).*(s2*f02).*real(fx.*conj(fy)) - 6*kappa*(sin( 4*psidpsi)*zrow).*(s4*f22).*real(fx.*conj(fy)) + 4*(cos(2*psidpsi)*zrow).*(s2*f02).*abs(fy).^2 + 5*kappa*(cos(2*psidpsi)*zrow).*(s2*f02).*abs(fy).^2 + 3*kappa*(cos(4*psidpsi)*zrow).*(s4*f22).*abs(fy).^2 + 4*(cos(2*psidpsi)*zrow).*(s2*f02).*abs(fz).^2 - 4*kappa*(cos( 2*psidpsi)*zrow).*(s2*f02).*abs(fz).^2);
    uy2 = -1/4*(cos(2*psidpsi)*zrow).*(s2*f02).*abs(fx).^2 - 1/8*kappa*(cos(2*psidpsi)*zrow).*(s2*f02).*abs(fx).^2 - 1/4*(cos(2*psidpsi)*zrow).*(s2*f11).*abs(fx).^2 - 1/8*kappa*(cos(2*psidpsi)*zrow).*(s2*f11).*abs(fx).^2 + 3/8*kappa*(cos(4*psidpsi)*zrow).*(s4*f22).*abs(fx).^2 + 3/4*kappa*(sin( 2*psidpsi)*zrow).*(s2*f02).*real(fx.*conj(fy)) + 3/4*kappa*(sin( 4*psidpsi)*zrow).*(s4*f22).*real(fx.*conj(fy)) - 1/4*(cos(2*psidpsi)*zrow).*(s2*f02).*abs(fy).^2 - 7/8*kappa*(cos(2*psidpsi)*zrow).*(s2*f02).*abs(fy).^2 - 1/4*(cos(2*psidpsi)*zrow).*(s2*f11).*abs(fy).^2 - 1/8*kappa*(cos(2*psidpsi)*zrow).*(s2*f11).*abs(fy).^2 - 3/8*kappa*(cos(4*psidpsi)*zrow).*(s4*f22).*abs(fy).^2 - 3/4*kappa*(cos(psidpsi)*zrow).*(s1*f12).*real(fx.*conj(fz)) + 3/4*kappa*(cos( 3*psidpsi)*zrow).*(s3*f12).*real(fx.*conj(fz)) + 3/2*kappa*(sin(psidpsi)*zrow).*(s1*f01).*real(fy.*conj(fz)) - 3/4*kappa*(sin(psidpsi)*zrow).*(s1*f12).*real(fy.*conj(fz)) + 3/4*kappa*(sin( 3*psidpsi)*zrow).*(s3*f12).*real(fy.*conj(fz)) - 1/4*(cos(2*psidpsi)*zrow).*(s2*f02).*abs(fz).^2 + kappa*(cos(2*psidpsi)*zrow).*(s2*f02).*abs(fz).^2 - 1/4*(cos(2*psidpsi)*zrow).*(s2*f11).*abs(fz).^2 + 1/4*kappa*(cos(2*psidpsi)*zrow).*(s2*f11).*abs(fz).^2 + (-1/4*(s0*f00).*abs(fx).^2 - 1/8*kappa*(s0*f00).*abs(fx).^2 + 1/4*(s0*f11).*abs(fx).^2 + 1/8*kappa*(s0*f11).*abs(fx).^2 - 1/4*(s0*f22).*abs(fx).^2 - 1/2*kappa*(s0*f22).*abs(fx).^2 - 1/4*(s0*f00).*abs(fy).^2 - 7/8*kappa*(s0*f00).*abs(fy).^2 + 1/4*(s0*f11).*abs(fy).^2 + 1/8*kappa*(s0*f11).*abs(fy).^2 - 1/4*(s0*f22).*abs(fy).^2 - 1/2*kappa*(s0*f22).*abs(fy).^2 - 1/4*(s0*f00).*abs(fz).^2 + kappa*(s0*f00).*abs(fz).^2 + 1/4*(s0*f11).*abs(fz).^2 - 1/4*kappa*(s0*f11).*abs(fz).^2 - 1/4*(s0*f22).*abs(fz).^2 + kappa*(s0*f22).*abs(fz).^2);
    uy3 = (3/16*kappa*(s0*f00).*abs(fx).^2 - 3/8*kappa*(s0*f11).*abs(fx).^2 + 3/8*kappa*(s0*f22).*abs(fx).^2 + 9/16*kappa*(s0*f00).*abs(fy).^2 - 3/8*kappa*(s0*f11).*abs(fy).^2 + 3/8*kappa*(s0*f22).*abs(fy).^2 - 3/4*kappa*(s0*f00).*abs(fz).^2 + 3/4*kappa*(s0*f11).*abs(fz).^2 - 3/4*kappa*(s0*f22).*abs(fz).^2) + 3/16*(kappa*(cos(2*psidpsi)*zrow).*(s2*f02).*abs(fx).^2 + 2*kappa*(cos(2*psidpsi)*zrow).*(s2*f11).*abs(fx).^2 - kappa*(cos(4*psidpsi)*zrow).*(s4*f22).*abs(fx).^2 - 2*kappa*(sin( 2*psidpsi)*zrow).*(s2*f02).*real(fx.*conj(fy)) - 2*kappa*(sin( 4*psidpsi)*zrow).*(s4*f22).*real(fx.*conj(fy)) + 3*kappa*(cos(2*psidpsi)*zrow).*(s2*f02).*abs(fy).^2 + 2*kappa*(cos(2*psidpsi)*zrow).*(s2*f11).*abs(fy).^2 + kappa*(cos(4*psidpsi)*zrow).*(s4*f22).*abs(fy).^2 + 4*kappa*(cos(psidpsi)*zrow).*(s1*f12).*real(fx.*conj(fz)) - 4*kappa*(cos(3*psidpsi)*zrow).*(s3*f12).*real(fx.*conj(fz)) - 8*kappa*(sin(psidpsi)*zrow).*(s1*f01).*real(fy.*conj(fz)) + 4*kappa*(sin(psidpsi)*zrow).*(s1*f12).*real(fy.*conj(fz)) - 4*kappa*(sin( 3*psidpsi)*zrow).*(s3*f12).*real(fy.*conj(fz)) - 4*kappa*(cos(2*psidpsi)*zrow).*(s2*f02).*abs(fz).^2 - 4*kappa*(cos(2*psidpsi)*zrow).*(s2*f11).*abs(fz).^2);

    % second focus
    fx = fxc2(:,:,1); fy = fyc2(:,:,1); fz = fzc2(:,:,1);
    for k=1:maxm
        fx = fx + fxc2(:,:,k+1)*cos(k*psi) + fxs2(:,:,k)*sin(k*psi);
        fy = fy + fyc2(:,:,k+1)*cos(k*psi) + fys2(:,:,k)*sin(k*psi);
        fz = fz + fzc2(:,:,k+1)*cos(k*psi) + fzs2(:,:,k)*sin(k*psi);
    end

    if sat2>0
        int = sat2*(abs(fx).^2+abs(fy).^2+abs(fz).^2);
        tmp = PulsedExcitation(int,1,pulse(1),pulse(2));
        tmp0 = tmp0 + triplet*tmp;
        tmp = sqrt(tmp./int);
        fx = fx.*tmp;
        fy = fy.*tmp;
        fz = fz.*tmp;
    else
        tmp0 = tmp0 + triplet*(abs(fx).^2+abs(fy).^2+abs(fz).^2);
    end
    
    tmp = (ux1*lam1 + ux2*lam2 + ux3*lam3)./tmp0/(4*maxm+1);
    volx1(:,:,1) = volx1(:,:,1) + tmp;
    for k = 1:maxm
        volx1(:,:,k+1) = volx1(:,:,k+1) + 2*tmp*cos(k*psi);
        volx1(:,:,maxm+1+k) = volx1(:,:,maxm+k+1) + 2*tmp*sin(k*psi);
    end

    tmp = (uy1*lam1 + uy2*lam2 + uy3*lam3)./tmp0/(4*maxm+1);
    voly1(:,:,1) = voly1(:,:,1) + tmp;
    for k = 1:maxm
        voly1(:,:,k+1) = voly1(:,:,k+1) + 2*tmp*cos(k*psi);
        voly1(:,:,maxm+1+k) = voly1(:,:,maxm+k+1) + 2*tmp*sin(k*psi);
    end

    tmp = sqrt(uv.^2+4*b.^2+4*b*uv*cos(psi));
    phi = real(acos((av^2-col*rhov.^2-tmp.^2*row)./(2*tmp*rhov)));
    s0 = 2*(pi-phi).*mat;
    s1 = 2*sin(phi).*mat;
    s2 = sin(2*phi).*mat;
    s3 = 2/3*sin(3*phi).*mat;
    s4 = 1/2*sin(4*phi).*mat;

    psidpsi = psi + real(acos((uv+2*b*cos(psi))./sqrt(uv.^2+4*b^2+4*b*uv*cos(psi))));
    psidpsi(~isfinite(psidpsi)) = 0;
    ux1 = (1/4*(s0*f00).*abs(fx).^2 + 5/16*kappa*(s0*f00).*abs(fx).^2 + 1/4*(s0*f22).*abs(fx).^2 + 1/8*kappa*(s0*f22).*abs(fx).^2 + 1/4*(s0*f00).*abs(fy).^2 - 1/16*kappa*(s0*f00).*abs(fy).^2 + 1/4*(s0*f22).*abs(fy).^2 + 1/8*kappa*(s0*f22).*abs(fy).^2 + 1/4*(s0*f00).*abs(fz).^2 - 1/4*kappa*(s0*f00).*abs(fz).^2 + 1/4*(s0*f22).*abs(fz).^2 - 1/4*kappa*(s0*f22).*abs(fz).^2) + 1/16*(-4*(cos(2*psidpsi)*zrow).*(s2*f02).*abs(fx).^2 - 5*kappa*(cos(2*psidpsi)*zrow).*(s2*f02).*abs(fx).^2 + 3*kappa*(cos(4*psidpsi)*zrow).*(s4*f22).*abs(fx).^2 - 6*kappa*(sin(2*psidpsi)*zrow).*(s2*f02).*real(fx.*conj(fy)) + 6*kappa*(sin(4*psidpsi)*zrow).*(s4*f22).*real(fx.*conj(fy)) - 4*(cos(2*psidpsi)*zrow).*(s2*f02).*abs(fy).^2 + kappa*(cos(2*psidpsi)*zrow).*(s2*f02).*abs(fy).^2 - 3*kappa*(cos(4*psidpsi)*zrow).*(s4*f22).*abs(fy).^2 - 4*(cos(2*psidpsi)*zrow).*(s2*f02).*abs(fz).^2 + 4*kappa*(cos(2*psidpsi)*zrow).*(s2*f02).*abs(fz).^2);
    ux2 = 1/4*(cos(2*psidpsi)*zrow).*(s2*f02).*abs(fx).^2 + 7/8*kappa*(cos(2*psidpsi)*zrow).*(s2*f02).*abs(fx).^2 + 1/4*(cos(2*psidpsi)*zrow).*(s2*f11).*abs(fx).^2 + 1/8*kappa*(cos(2*psidpsi)*zrow).*(s2*f11).*abs(fx).^2 - 3/8*kappa*(cos(4*psidpsi)*zrow).*(s4*f22).*abs(fx).^2 + 3/4*kappa*(sin(2*psidpsi)*zrow).*(s2*f02).*real(fx.*conj(fy)) - 3/4*kappa*(sin(4*psidpsi)*zrow).*(s4*f22).*real(fx.*conj(fy)) + 1/4*(cos(2*psidpsi)*zrow).*(s2*f02).*abs(fy).^2 + 1/8*kappa*(cos(2*psidpsi)*zrow).*(s2*f02).*abs(fy).^2 + 1/4*(cos(2*psidpsi)*zrow).*(s2*f11).*abs(fy).^2 + 1/8*kappa*(cos(2*psidpsi)*zrow).*(s2*f11).*abs(fy).^2 + 3/8*kappa*(cos(4*psidpsi)*zrow).*(s4*f22).*abs(fy).^2 + 3/2*kappa*(cos(psidpsi)*zrow).*(s1*f01).*real(fx.*conj(fz)) - 3/4*kappa*(cos(psidpsi)*zrow).*(s1*f12).*real(fx.*conj(fz)) - 3/4*kappa*(cos(3*psidpsi)*zrow).*(s3*f12).*real(fx.*conj(fz)) - 3/4*kappa*(sin(psidpsi)*zrow).*(s1*f12).*real(fy.*conj(fz)) - 3/4*kappa*(sin(3*psidpsi)*zrow).*(s3*f12).*real(fy.*conj(fz)) + 1/4*(cos(2*psidpsi)*zrow).*(s2*f02).*abs(fz).^2 - kappa*(cos(2*psidpsi)*zrow).*(s2*f02).*abs(fz).^2 + 1/4*(cos(2*psidpsi)*zrow).*(s2*f11).*abs(fz).^2 - 1/4*kappa*(cos(2*psidpsi)*zrow).*(s2*f11).*abs(fz).^2 + (-1/4*(s0*f00).*abs(fx).^2 - 7/8*kappa*(s0*f00).*abs(fx).^2 + 1/4*(s0*f11).*abs(fx).^2 + 1/8*kappa*(s0*f11).*abs(fx).^2 - 1/4*(s0*f22).*abs(fx).^2 - 1/2*kappa*(s0*f22).*abs(fx).^2 - 1/4*(s0*f00).*abs(fy).^2 - 1/8*kappa*(s0*f00).*abs(fy).^2 + 1/4*(s0*f11).*abs(fy).^2 + 1/8*kappa*(s0*f11).*abs(fy).^2 - 1/4*(s0*f22).*abs(fy).^2 - 1/2*kappa*(s0*f22).*abs(fy).^2 - 1/4*(s0*f00).*abs(fz).^2 + kappa*(s0*f00).*abs(fz).^2 + 1/4*(s0*f11).*abs(fz).^2 - 1/4*kappa*(s0*f11).*abs(fz).^2 - 1/4*(s0*f22).*abs(fz).^2 + kappa*(s0*f22).*abs(fz).^2);
    ux3 = (9/16*kappa*(s0*f00).*abs(fx).^2 - 3/8*kappa*(s0*f11).*abs(fx).^2 + 3/8*kappa*(s0*f22).*abs(fx).^2 + 3/16*kappa*(s0*f00).*abs(fy).^2 - 3/8*kappa*(s0*f11).*abs(fy).^2 + 3/8*kappa*(s0*f22).*abs(fy).^2 - 3/4*kappa*(s0*f00).*abs(fz).^2 + 3/4*kappa*(s0*f11).*abs(fz).^2 - 3/4*kappa*(s0*f22).*abs(fz).^2) - 3/16*(3*kappa*(cos(2*psidpsi)*zrow).*(s2*f02).*abs(fx).^2 + 2*kappa*(cos(2*psidpsi)*zrow).*(s2*f11).*abs(fx).^2 - kappa*(cos(4*psidpsi)*zrow).*(s4*f22).*abs(fx).^2 + 2*kappa*(sin(2*psidpsi)*zrow).*(s2*f02).*real(fx.*conj(fy)) - 2*kappa*(sin(4*psidpsi)*zrow).*(s4*f22).*real(fx.*conj(fy)) + kappa*(cos(2*psidpsi)*zrow).*(s2*f02).*abs(fy).^2 + 2*kappa*(cos(2*psidpsi)*zrow).*(s2*f11).*abs(fy).^2 + kappa*(cos(4*psidpsi)*zrow).*(s4*f22).*abs(fy).^2 + 8*kappa*(cos(psidpsi)*zrow).*(s1*f01).*real(fx.*conj(fz)) - 4*kappa*(cos(psidpsi)*zrow).*(s1*f12).*real(fx.*conj(fz)) - 4*kappa*(cos(3*psidpsi)*zrow).*(s3*f12).*real(fx.*conj(fz)) - 4*kappa*(sin(psidpsi)*zrow).*(s1*f12).*real(fy.*conj(fz)) - 4*kappa*(sin(3*psidpsi)*zrow).*(s3*f12).*real(fy.*conj(fz)) - 4*kappa*(cos(2*psidpsi)*zrow).*(s2*f02).*abs(fz).^2 - 4*kappa*(cos(2*psidpsi)*zrow).*(s2*f11).*abs(fz).^2);

    tmp = (ux1*lam1 + ux2*lam2 + ux3*lam3)./tmp0/(4*maxm+1);
    volx2(:,:,1) = volx2(:,:,1) + tmp;
    for k = 1:maxm
        volx2(:,:,k+1) = volx2(:,:,k+1) + 2*tmp*cos(k*psi);
        volx2(:,:,maxm+1+k) = volx2(:,:,maxm+k+1) + 2*tmp*sin(k*psi);
    end

    tmp = sqrt(uv.^2+b.^2+2*b*uv*cos(psi));
    phi = real(acos((av^2-col*rhov.^2-tmp.^2*row)./(2*tmp*rhov)));
    s0 = 2*(pi-phi).*mat;
    s1 = 2*sin(phi).*mat;
    s2 = sin(2*phi).*mat;
    s3 = 2/3*sin(3*phi).*mat;
    s4 = 1/2*sin(4*phi).*mat;

    psidpsi = psi + real(acos((uv+b*cos(psi))./sqrt(uv.^2+b^2+2*b*uv*cos(psi))));
    psidpsi(~isfinite(psidpsi)) = 0;
    uy1 = (1/4*(s0*f00).*abs(fx).^2 - 1/16*kappa*(s0*f00).*abs(fx).^2 + 1/4*(s0*f22).*abs(fx).^2 + 1/8*kappa*(s0*f22).*abs(fx).^2 + 1/4*(s0*f00).*abs(fy).^2 + 5/16*kappa*(s0*f00).*abs(fy).^2 + 1/4*(s0*f22).*abs(fy).^2 + 1/8*kappa*(s0*f22).*abs(fy).^2 + 1/4*(s0*f00).*abs(fz).^2 - 1/4*kappa*(s0*f00).*abs(fz).^2 + 1/4*(s0*f22).*abs(fz).^2 - 1/4*kappa*(s0*f22).*abs(fz).^2) + 1/16*(4*(cos(2*psidpsi)*zrow).*(s2*f02).*abs(fx).^2 - kappa*(cos(2*psidpsi)*zrow).*(s2*f02).*abs(fx).^2 - 3*kappa*(cos(4*psidpsi)*zrow).*(s4*f22).*abs(fx).^2 - 6*kappa*(sin( 2*psidpsi)*zrow).*(s2*f02).*real(fx.*conj(fy)) - 6*kappa*(sin( 4*psidpsi)*zrow).*(s4*f22).*real(fx.*conj(fy)) + 4*(cos(2*psidpsi)*zrow).*(s2*f02).*abs(fy).^2 + 5*kappa*(cos(2*psidpsi)*zrow).*(s2*f02).*abs(fy).^2 + 3*kappa*(cos(4*psidpsi)*zrow).*(s4*f22).*abs(fy).^2 + 4*(cos(2*psidpsi)*zrow).*(s2*f02).*abs(fz).^2 - 4*kappa*(cos( 2*psidpsi)*zrow).*(s2*f02).*abs(fz).^2);
    uy2 = -1/4*(cos(2*psidpsi)*zrow).*(s2*f02).*abs(fx).^2 - 1/8*kappa*(cos(2*psidpsi)*zrow).*(s2*f02).*abs(fx).^2 - 1/4*(cos(2*psidpsi)*zrow).*(s2*f11).*abs(fx).^2 - 1/8*kappa*(cos(2*psidpsi)*zrow).*(s2*f11).*abs(fx).^2 + 3/8*kappa*(cos(4*psidpsi)*zrow).*(s4*f22).*abs(fx).^2 + 3/4*kappa*(sin( 2*psidpsi)*zrow).*(s2*f02).*real(fx.*conj(fy)) + 3/4*kappa*(sin( 4*psidpsi)*zrow).*(s4*f22).*real(fx.*conj(fy)) - 1/4*(cos(2*psidpsi)*zrow).*(s2*f02).*abs(fy).^2 - 7/8*kappa*(cos(2*psidpsi)*zrow).*(s2*f02).*abs(fy).^2 - 1/4*(cos(2*psidpsi)*zrow).*(s2*f11).*abs(fy).^2 - 1/8*kappa*(cos(2*psidpsi)*zrow).*(s2*f11).*abs(fy).^2 - 3/8*kappa*(cos(4*psidpsi)*zrow).*(s4*f22).*abs(fy).^2 - 3/4*kappa*(cos(psidpsi)*zrow).*(s1*f12).*real(fx.*conj(fz)) + 3/4*kappa*(cos( 3*psidpsi)*zrow).*(s3*f12).*real(fx.*conj(fz)) + 3/2*kappa*(sin(psidpsi)*zrow).*(s1*f01).*real(fy.*conj(fz)) - 3/4*kappa*(sin(psidpsi)*zrow).*(s1*f12).*real(fy.*conj(fz)) + 3/4*kappa*(sin( 3*psidpsi)*zrow).*(s3*f12).*real(fy.*conj(fz)) - 1/4*(cos(2*psidpsi)*zrow).*(s2*f02).*abs(fz).^2 + kappa*(cos(2*psidpsi)*zrow).*(s2*f02).*abs(fz).^2 - 1/4*(cos(2*psidpsi)*zrow).*(s2*f11).*abs(fz).^2 + 1/4*kappa*(cos(2*psidpsi)*zrow).*(s2*f11).*abs(fz).^2 + (-1/4*(s0*f00).*abs(fx).^2 - 1/8*kappa*(s0*f00).*abs(fx).^2 + 1/4*(s0*f11).*abs(fx).^2 + 1/8*kappa*(s0*f11).*abs(fx).^2 - 1/4*(s0*f22).*abs(fx).^2 - 1/2*kappa*(s0*f22).*abs(fx).^2 - 1/4*(s0*f00).*abs(fy).^2 - 7/8*kappa*(s0*f00).*abs(fy).^2 + 1/4*(s0*f11).*abs(fy).^2 + 1/8*kappa*(s0*f11).*abs(fy).^2 - 1/4*(s0*f22).*abs(fy).^2 - 1/2*kappa*(s0*f22).*abs(fy).^2 - 1/4*(s0*f00).*abs(fz).^2 + kappa*(s0*f00).*abs(fz).^2 + 1/4*(s0*f11).*abs(fz).^2 - 1/4*kappa*(s0*f11).*abs(fz).^2 - 1/4*(s0*f22).*abs(fz).^2 + kappa*(s0*f22).*abs(fz).^2);
    uy3 = (3/16*kappa*(s0*f00).*abs(fx).^2 - 3/8*kappa*(s0*f11).*abs(fx).^2 + 3/8*kappa*(s0*f22).*abs(fx).^2 + 9/16*kappa*(s0*f00).*abs(fy).^2 - 3/8*kappa*(s0*f11).*abs(fy).^2 + 3/8*kappa*(s0*f22).*abs(fy).^2 - 3/4*kappa*(s0*f00).*abs(fz).^2 + 3/4*kappa*(s0*f11).*abs(fz).^2 - 3/4*kappa*(s0*f22).*abs(fz).^2) + 3/16*(kappa*(cos(2*psidpsi)*zrow).*(s2*f02).*abs(fx).^2 + 2*kappa*(cos(2*psidpsi)*zrow).*(s2*f11).*abs(fx).^2 - kappa*(cos(4*psidpsi)*zrow).*(s4*f22).*abs(fx).^2 - 2*kappa*(sin( 2*psidpsi)*zrow).*(s2*f02).*real(fx.*conj(fy)) - 2*kappa*(sin( 4*psidpsi)*zrow).*(s4*f22).*real(fx.*conj(fy)) + 3*kappa*(cos(2*psidpsi)*zrow).*(s2*f02).*abs(fy).^2 + 2*kappa*(cos(2*psidpsi)*zrow).*(s2*f11).*abs(fy).^2 + kappa*(cos(4*psidpsi)*zrow).*(s4*f22).*abs(fy).^2 + 4*kappa*(cos(psidpsi)*zrow).*(s1*f12).*real(fx.*conj(fz)) - 4*kappa*(cos(3*psidpsi)*zrow).*(s3*f12).*real(fx.*conj(fz)) - 8*kappa*(sin(psidpsi)*zrow).*(s1*f01).*real(fy.*conj(fz)) + 4*kappa*(sin(psidpsi)*zrow).*(s1*f12).*real(fy.*conj(fz)) - 4*kappa*(sin( 3*psidpsi)*zrow).*(s3*f12).*real(fy.*conj(fz)) - 4*kappa*(cos(2*psidpsi)*zrow).*(s2*f02).*abs(fz).^2 - 4*kappa*(cos(2*psidpsi)*zrow).*(s2*f11).*abs(fz).^2);

    tmp = (uy1*lam1 + uy2*lam2 + uy3*lam3)./tmp0/(4*maxm+1);
    voly2(:,:,1) = voly2(:,:,1) + tmp;
    for k = 1:maxm
        voly2(:,:,k+1) = voly2(:,:,k+1) + 2*tmp*cos(k*psi);
        voly2(:,:,maxm+1+k) = voly2(:,:,maxm+k+1) + 2*tmp*sin(k*psi);
    end

end
warning on MATLAB:divideByZero

mdf.NA = NA;
mdf.n0 = n0;
mdf.n = n;
mdf.n1 = n1;
mdf.focpos = focpos;
mdf.lambda = lamem;
mdf.mag = mag;
mdf.av = av;
mdf.zpin = zpin;
mdf.atf = atf;
mdf.kappa = kappa;
mdf.lt = lt;
mdf.pulse = pulse;
mdf.sat = sat;
mdf.triplet = triplet;
mdf.volx1 = volx1;
mdf.voly1 = voly1;
mdf.volx2 = volx2;
mdf.voly2 = voly2;
mdf.qv = qv;
mdf.qp = qp;



