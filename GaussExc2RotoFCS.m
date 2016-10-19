function [ux, uy, pp, po, op, oo] = GaussExc2RotoFCS(exc, lamem, mag, av, zpin, atf, gamma)

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
[m0, m1, m2, m3, m4] = RotoDiffCoef(gamma); % Calculated from CEFRotoFCS.nb
mm = 0*m4;
for L=0:4
    eval(['mm(1:L+1,1:2*L+1,:,:,L+1) = m' int2str(L) ';'])
end

poyntx = zeros(size(f00,1),size(f00,2),5,5,9);
poynty = poyntx;

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
    
    ss = {'xx','yy','zz','xy','xz','yz'};
    
    for L=0:4
        for M=0:L
            for K=-L:L
                for jh=1:6
                    for jf=1:6
                        poyntx(:,:,L+1,M+1,L+K+1) = poyntx(:,:,L+1,M+1,L+K+1) + mm(M+1,L+K+1,jh,jf,L+1)*h(:,:,jh).*f(:,:,jf);
                        poynty(:,:,L+1,M+1,L+K+1) = poynty(:,:,L+1,M+1,L+K+1) + mm(M+1,L+K+1,jh,jf,L+1)*h(:,:,jh).*g(:,:,jf);                        
                    end
                end
            end
        end
    end
        
    ux = ux + real(sum(rho(:,1)'*poyntx(:,:,1,1,1)));
    uy = uy + real(sum(rho(:,1)'*poynty(:,:,1,1,1)));
    for L = 0:4
        for M = 0:L
            for K = -L:L
                pp(L+1,M+1) = pp(L+1,M+1) + (1+sign(M))*sum(rho(:,1)'*real(poyntx(:,:,L+1,M+1,L+K+1).*conj(poyntx(:,:,L+1,M+1,L+K+1))));
                po(L+1,M+1) = po(L+1,M+1) + (1+sign(M))*sum(rho(:,1)'*real(poyntx(:,:,L+1,M+1,L+K+1).*conj(poynty(:,:,L+1,M+1,L+K+1))));
                op(L+1,M+1) = op(L+1,M+1) + (1+sign(M))*sum(rho(:,1)'*real(poynty(:,:,L+1,M+1,L+K+1).*conj(poynty(:,:,L+1,M+1,L+K+1))));
                oo(L+1,M+1) = oo(L+1,M+1) + (1+sign(M))*sum(rho(:,1)'*real(1i^M*poyntx(:,:,L+1,M+1,L+K+1).*conj(poyntx(:,:,L+1,M+1,L+K+1))));
            end
        end
    end
    jpsi
end

return

% the following is no longer valid!
for L=0:2:4
    tmpp(L/2+1) = pp(L+1,1) + sum(pp(L+1,2:end));
    if 0 % Anastasia mode
        tmpo(L/2+1) = po(L+1,1) + sum(real(1i.^(1:6).*po(L+1,2:end)));
        tmop(L/2+1) = op(L+1,1) + sum(real(1i.^(1:6).*op(L+1,2:end)));
    else % Christoph mode
        tmpo(L/2+1) = po(L+1,1) + sum(real(po(L+1,2:end)));
        tmop(L/2+1) = op(L+1,1) + sum(real(op(L+1,2:end)));
    end
end
tmpp = tmpp/po(1);
tmpo = tmpo/po(1);
tmop = tmop/po(1);

drot = 1/20/6;
t = 0:100;
plot(t,tmpp(1)+tmpp(2)*exp(-6*drot*t)+tmpp(3)*exp(-20*drot*t),...
    t,tmpo(1)+tmpo(2)*exp(-6*drot*t)+tmpo(3)*exp(-20*drot*t),...
    t,tmop(1)+tmop(2)*exp(-6*drot*t)+tmop(3)*exp(-20*drot*t))
xlabel('time (ns)');
ylabel('correlation (a.u.)')

return

drot = 1/20/6;
d0 = 0;
d1 = (1/20-1/40)/6;
t = 0:100;
zpp0 = 0*t;
zpp1 = zpp0;
zpo0 = zpp0;
zpo1 = zpp0;
zop0 = zpp0;
zop1 = zpp0;
for L=0:6
    for M=-L:L
        zpp0 = zpp0 + real(pp(L+1,abs(M)+1))*exp(-L*(L+1)*drot*t - M^2*d0*t);
        zpo0 = zpo0 + real(po(L+1,abs(M)+1))*exp(-L*(L+1)*drot*t - M^2*d0*t);
        zop0 = zop0 + real(op(L+1,abs(M)+1))*exp(-L*(L+1)*drot*t - M^2*d0*t);        

        zpp1 = zpp1 + real(pp(L+1,abs(M)+1))*exp(-L*(L+1)*drot*t - M^2*d1*t);
        zpo1 = zpo1 + real(po(L+1,abs(M)+1))*exp(-L*(L+1)*drot*t - M^2*d1*t);
        zop1 = zop1 + real(op(L+1,abs(M)+1))*exp(-L*(L+1)*drot*t - M^2*d1*t);        
    end
end
plot(t,zpp0,t,zpp1,t,zpo0,t,zpo1,t,zop0,t,zop1)
xlabel('time (ns)');
ylabel('correlation (a.u.)')

