function [exci, emix, emiy] = GaussExc2Aniso(exc, lamem, mag, av, zpin, atf, gamma)

ni = 1.0; % it is assumed that imaging medium is air

NA = exc.NA;
n0 = exc.n0;
n = exc.n;
n1 = exc.n1;
focpos = exc.focpos(1);
maxm = exc.maxm;
rho = exc.rho;
z = exc.z;
fxc = exc.fxc;
fxs = exc.fxs;
fyc = exc.fyc;
fys = exc.fys;
fzc = exc.fzc;
fzs = exc.fzs;
maxnum = exc.maxnum;

ki = 2*pi/lamem*ni;
k0 = 2*pi/lamem*n0; 

if length(rho)>1
    drho = rho(2,1)-rho(1,1);
    rhov = mag*(rho(1,1):drho:(rho(end,1)+av/mag));
    uv = mag*rho(:,1);
else
    if length(z)>1
        drho = z(2)-z(1);
        rhov = mag*(rho:drho:(rho+av/mag));
        uv = mag*rho;
    else
        rho = 1;
        rhov = 0;
        uv = 0;
    end
end
zv = z(1,:);

chimin = 0;
chimax = asin(NA/ni/mag); % excluding emission above the critical angle
dchi = chimax/maxnum;
chi = (chimin:dchi:chimax)';
if isempty(chi) 
    chi=0; 
end
ci = cos(chi);
si = sin(chi);
s0 = ni/n0(1)*mag*si; % Abbe's sine condition
c0 = sqrt(1-s0.^2);

d0 = 2*pi/lamem*exc.d0; d = 2*pi/lamem*exc.d; d1 = 2*pi/lamem*exc.d1;
[v,pc,ps] = DipoleL(asin(s0),2*pi*zv/lamem,n0,n,n1,d0,d,d1);
row = ones(size(zv));

if nargin>5 && ~isempty(atf)
    if length(atf)==1 % accounting for reflection losses when using water immersion; atf = ref index of cover slide
        [~, ~, tmp, tms] = Fresnel(n0(1)*c0,n0(1),atf);
        [~, ~, mp, ms] = Fresnel(sqrt(atf^2-n0(1)^2*s0.^2),atf,n0(1));
    else % + aberration for water immersion; atf(2) = thickness mismatch
        [~, ~, tmp, tms] = Fresnel(n0(1)*c0,[n0(1) atf(1) atf(1)], 2*pi*atf(2)/lamem);
        [~, ~, mp, ms] = Fresnel(sqrt(atf(1)^2-n0(1)^2*s0.^2),[atf(1) n0(1) n0(1)], -2*pi*atf(2)/lamem);
    end
    v = ((tmp.*mp).'*row).*v;
    pc = ((tmp.*mp).'*row).*pc;
    ps = ((tms.*ms).'*row).*ps;
end

phase = rem(-k0(1)*c0*focpos + ki*zpin*ci,2*pi);
ez = exp(1i*phase);

if length(chi)>1
    fac = dchi*si.*sqrt(ci./c0);
else
    fac = 1;
end
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

col = ones(size(uv)); 
row = ones(size(rhov)); 
mat = ones(size(uv))*rhov; 
if mat==0
    mat = 1;
end
phi = real(acos((av^2-col*rhov.^2-uv.^2*row)./(2*uv*rhov)));
a0 = 2*(pi-phi).*mat;
a1 = 2*sin(phi).*mat;
a2 = sin(2*phi).*mat;
a3 = 2/3*sin(3*phi).*mat;
a4 = 1/2*sin(4*phi).*mat;
f00 = a0*f00;
f01 = a1*f01;
f02 = a2*f02;
f112 = a2*f11;
f11 = a0*f11;
f224 = a4*f22;
f22 = a0*f22;
f123 = a3*f12;
f12 = a1*f12;

ux = zeros(1,1);
uy = ux; 
pp = zeros(4+1,4+1); 
po = pp; 
op = pp;
oo = pp;

if nargin<7 || isempty(gamma)
    gamma = 0;
end
[m0ex, m2ex] = RotoAnisoCoef; 
[m0em, m2em] = RotoAnisoCoef(gamma);

exci = zeros(size(f00,1),size(f00,2),2,5,9);
emix = exci;
emiy = exci;

for jpsi=0:2*maxm
    psi = jpsi/(2*maxm+1)*2*pi;
    
    hx = fxc(:,:,1); hy = fyc(:,:,1); hz = fzc(:,:,1);
    for k=1:maxm
        hx = hx + fxc(:,:,k+1)*cos(k*psi) + fxs(:,:,k)*sin(k*psi);
        hy = hy + fyc(:,:,k+1)*cos(k*psi) + fys(:,:,k)*sin(k*psi);
        hz = hz + fzc(:,:,k+1)*cos(k*psi) + fzs(:,:,k)*sin(k*psi);
    end
    
    h(:,:,1) = abs(hx).^2; 
    h(:,:,2) = abs(hy).^2; 
    h(:,:,3) = abs(hz).^2; 
    h(:,:,4) = 2*real(hx.*conj(hy)); 
    h(:,:,5) = 2*real(hx.*conj(hz));
    h(:,:,6) = 2*real(hy.*conj(hz));
    f(:,:,1) = f00 + 0.5*f22 - f02*cos(2*psi) + 0.5*f224*cos(4*psi);
    f(:,:,2) = 0.5*(f22 - f224*cos(4*psi));
    f(:,:,3) = 0.5*(f11 + f112*cos(2*psi));
    f(:,:,4) = f02*sin(2*psi) - f224*sin(4*psi);
    f(:,:,5) = (f01 - 0.5*f12)*cos(psi) - 0.5*f123*cos(3*psi);
    f(:,:,6) = 0.5*(f12*sin(psi) + f123*sin(3*psi));
    psi = psi - pi/2;
    g(:,:,2) = f00 + 0.5*f22 - f02*cos(2*psi) + 0.5*f224*cos(4*psi);
    g(:,:,1) = 0.5*(f22 - f224*cos(4*psi));
    g(:,:,3) = 0.5*(f11 + f112*cos(2*psi));
    g(:,:,4) = - f02*sin(2*psi) + f224*sin(4*psi);
    g(:,:,6) = (f01 - 0.5*f12)*cos(psi) - 0.5*f123*cos(3*psi);
    g(:,:,5) = - 0.5*(f12*sin(psi) + f123*sin(3*psi));
    
    % ss = {'xx','yy','zz','xy','xz','yz'};

    for j=1:6
        exci(:,:,1,1,1) = exci(:,:,1,1,1) + m0ex(j)*h(:,:,j);
        emix(:,:,1,1,1) = emix(:,:,1,1,1) + m0em(j)*f(:,:,j);
        emiy(:,:,1,1,1) = emiy(:,:,1,1,1) + m0em(j)*g(:,:,j);
    end
    
    for M=0:2
        for K=-2:2
            for j=1:6
                exci(:,:,2,M+1,2+K+1) = exci(:,:,2,M+1,2+K+1) + m2ex(M+1,2+K+1,j)*h(:,:,j);
                emix(:,:,2,M+1,2+K+1) = emix(:,:,2,M+1,2+K+1) + m2em(M+1,2+K+1,j)*f(:,:,j);
                emiy(:,:,2,M+1,2+K+1) = emiy(:,:,2,M+1,2+K+1) + m2em(M+1,2+K+1,j)*g(:,:,j);
            end
        end
    end
        
    jpsi
end

