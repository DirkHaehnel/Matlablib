function [ux, uy, pp, po, op] = GaussExc2RotoFCS(exc, lamem, mag, av, zpin, atf, sat)

maxL = 6;
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

ux = zeros(maxL+1,maxL+1);
uy = ux; 
pp = ux; 
po = ux; 
op = ux;
tmpxx = zeros(size(rho,1),size(z,2),maxL+1,maxL+1); 
tmpxy = tmpxx;
if ~(nargin<12 || isempty(sat) || sat==0)
    sat = sat/max(abs(fxc(1,:,1)).^2 + abs(fyc(1,:,1)).^2 + abs(fzc(1,:,1)).^2);
end

warning off MATLAB:divideByZero
[weight,x] = gausslegendrecof(2*maxL);
alv = (0.5:2*maxL)/(2*maxL)*2*pi; bev = acos(x); weight = weight*mean(diff(alv));
%clear tst; for j=0:maxL for k=0:j tst(j+1,k+1)=sum(weight'*(SphericalHarmonic(j,k,bev,alv).*conj(SphericalHarmonic(j,k,bev,alv)))); end; end
for j=0:2*maxm
    psi = j/(2*maxm+1)*2*pi;

    hx = fxc(:,:,1); hy = fyc(:,:,1); hz = fzc(:,:,1);
    for k=1:maxm
        hx = hx + fxc(:,:,k+1)*cos(k*psi) + fxs(:,:,k)*sin(k*psi);
        hy = hy + fyc(:,:,k+1)*cos(k*psi) + fys(:,:,k)*sin(k*psi);
        hz = hz + fzc(:,:,k+1)*cos(k*psi) + fzs(:,:,k)*sin(k*psi);
    end
    
    hxx = abs(hx).^2;
    hyy = abs(hy).^2;
    hzz = abs(hz).^2;
    hxy = 2*real(hx*conj(hy));
    hxz = 2*real(hx*conj(hz));
    hyz = 2*real(hy*conj(hz));
    fxx = f00 + 0.5*f22 - f02*cos(2*psi) + 0.5*f224*cos(4*psi);
    fyy = 0.5*(f22 - ff24*cos(4*psi));
    fzz = 0.5*(f11 + f112*cos(2*psi));
    fxy = f02*sin(2*psi) - f224*sin(4*psi);
    fxz = (f01 - 0.5*f12)*cos(psi) - 0.5*f123*cos(3*psi);
    fyz = 0.5*(f12*sin(psi) + f123*sin(3*psi));
    psi = psi - pi/2;
    gxx = f00 + 0.5*f22 - f02*cos(2*psi) + 0.5*f224*cos(4*psi);
    gyy = 0.5*(f22 - ff24*cos(4*psi));
    gzz = 0.5*(f11 + f112*cos(2*psi));
    gxy = f02*sin(2*psi) - f224*sin(4*psi);
    gxz = (f01 - 0.5*f12)*cos(psi) - 0.5*f123*cos(3*psi);
    gyz = 0.5*(f12*sin(psi) + f123*sin(3*psi));

    for jal=1:length(alv)
        al = alv(jal);
        for jbe=1:length(bev)
            be = bev(jbe);
            vx = sin(be)*cos(al);
            vy = sin(be)*sin(al);
            vz = cos(be);
            for j=0:4
                for m = -j:j
                    for k = -j:j
            for j1=0:4
                for m1 = -j1:j1
                    for k1 = -j1:j1
            for j2=0:4
                for m2 = -j2:j2
                    for k2 = -j2:j2
                        tmp = clebsch(j1,m1,j2,m2,j,-m)*clebsch(j1,k1,j2,k2,j,-k);
                        if tmp>0
                            htmp = WignerCoefficients(j2,m2,k2,hxx,hyy,hzz,hxy,hxz,hyz,vx,vy,vz);
                            poyntx(:,:,j+1,j+1+m,j+1+k) = poyntx(:,:,j+1,j+1+m,j+1+k) + ...
                                tmp*WignerCoefficients(j1,m1,k1,fxx,fyy,fzz,fxy,fxz,fyz,vx,vy,vz).*htmp;
                            poynty(:,:,j+1,j+1+m,j+1+k) = poynty(:,:,j+1,j+1+m,j+1+k) + ...
                                tmp*WignerCoefficients(j1,m1,k1,gxx,gyy,gzz,gxy,gxz,gyz,vy,-vx,vz).*htmp;
                        end
                    end
                end
            end
                    end
                end
            end
                    end
                end
            end
            
            for L = 0:maxL
                for M = 0:L
                    tmpxx(:,:,L+1,M+1) = tmpxx(:,:,L+1,M+1) + (tmp.*poyntx)*conj(SphericalHarmonic(L,M,be,al));
                    tmpxy(:,:,L+1,M+1) = tmpxy(:,:,L+1,M+1) + (tmp.*poynty)*conj(SphericalHarmonic(L,M,be,al));
                end
            end
            disp(['beta index = ' int2str(jbe)])
        end
        disp(['alpha index = ' int2str(jal)])
    end
    disp(['psi index = ' int2str(j)])
    for L = 0:maxL
        for M = 0:L
            ux(L+1,M+1) = ux(L+1,M+1) + sum(rho(:,1)'*tmpxx(:,:,L+1,M+1));
            uy(L+1,M+1) = uy(L+1,M+1) + sum(rho(:,1)'*tmpxy(:,:,L+1,M+1));
            pp(L+1,M+1) = pp(L+1,M+1) + sum(rho(:,1)'*real((tmpxx(:,:,L+1,M+1).*conj(tmpxy(:,:,L+1,M+1)))));
            po(L+1,M+1) = po(L+1,M+1) + sum(rho(:,1)'*abs(tmpxx(:,:,L+1,M+1)).^2);
            op(L+1,M+1) = op(L+1,M+1) + sum(rho(:,1)'*abs(tmpxy(:,:,L+1,M+1)).^2);
        end
    end
end

warning on MATLAB:divideByZero

return

for L=0:2:4
    tmpp(L/2+1) = pp(L+1,1) + 2*sum(pp(L+1,2:end));
    if 0 % Anastasia mode
        tmpo(L/2+1) = po(L+1,1) + 2*sum(real(1i.^(1:6).*po(L+1,2:end)));
        tmop(L/2+1) = op(L+1,1) + 2*sum(real(1i.^(1:6).*op(L+1,2:end)));
    else % Christoph mode
        tmpo(L/2+1) = po(L+1,1) + 2*sum(real(po(L+1,2:end)));
        tmop(L/2+1) = op(L+1,1) + 2*sum(real(op(L+1,2:end)));
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

% expansion of fxx*vx^2+fyy*vy^2+fzz*vz^2+fxy*vx*vy+fxz*vx*vz+fyz*vy*vz 
% in Wigner eigenfunctions of SO(3) group
function f = WignerCoefficients(j,m,k,fxx,fyy,fzz,fxy,fxz,fyz,vx,vy,vz)

switch j
    case 0
        f = (2.*sqrt(2).*(fxx + fyy + fzz).*pi.*(vx.^2 + vy.^2 + vz.^2))/3;

    case 1
        f = 0;
        
    case 2
        switch m
            case -2
                switch k
                    case -2
                        f = (fxx + 1i.*fxy - fyy).*pi.*(vx - 1i.*vy).^2/sqrt(10);
                    case -1
                        f = ((fxz + 1i.*fyz).*pi.*(vx - 1i.*vy).^2)/sqrt(10);
                    case 0
                        f = -(((fxx + fyy - 2.*fzz).*pi.*(vx - 1i.*vy).^2)/sqrt(15));
                    case 1
                        f = -(((fxz - 1i.*fyz).*pi.*(vx - 1i.*vy).^2)/sqrt(10));
                    case 2
                        f =  ((fxx - 1i.*fxy - fyy).*pi.*(vx - 1i.*vy).^2)/sqrt(10);
                end
            case -1
                switch k
                    case -2
                        f = sqrt(2/5).*(fxx + 1i.*fxy - fyy).*pi.*(vx - 1i.*vy).*vz;
                    case -1
                        f = sqrt(2/5).*(fxz + 1i.*fyz).*pi.*(vx - 1i.*vy).*vz;
                    case 0
                        f = (-2.*(fxx + fyy - 2.*fzz).*pi.*(vx - 1i.*vy).*vz)/sqrt(15);
                    case 1
                        f = -(sqrt(2/5).*(fxz - 1i.*fyz).*pi.*(vx - 1i.*vy).*vz);
                    case 2
                        f =  sqrt(2/5).*(fxx - 1i.*fxy - fyy).*pi.*(vx - 1i.*vy).*vz;
                end
            case 0
                switch k
                    case -2
                        f = -(((fxx + 1i.*fxy - fyy).*pi.*(vx.^2 + vy.^2 - 2.*vz.^2))/sqrt(15));
                    case -1
                        f = -(((fxz + 1i.*fyz).*pi.*(vx.^2 + vy.^2 - 2.*vz.^2))/sqrt(15));
                    case 0
                        f = (sqrt(2/5).*(fxx + fyy - 2.*fzz).*pi.*(vx.^2 + vy.^2 - 2.*vz.^2))/3;
                    case 1
                        f =  ((fxz - 1i.*fyz).*pi.*(vx.^2 + vy.^2 - 2.*vz.^2))/sqrt(15);
                    case 2
                        f = -(((fxx - 1i.*fxy - fyy).*pi.*(vx.^2 + vy.^2 - 2.*vz.^2))/sqrt(15));
                end
            case 1
                switch k
                    case -2
                        f = -(sqrt(2/5).*(fxx + 1i.*fxy - fyy).*pi.*(vx + 1i.*vy).*vz);
                    case -1
                        f = -(sqrt(2/5).*(fxz + 1i.*fyz).*pi.*(vx + 1i.*vy).*vz);
                    case 0
                        f =  (2.*(fxx + fyy - 2.*fzz).*pi.*(vx + 1i.*vy).*vz)/sqrt(15);
                    case 1
                        f = sqrt(2/5).*(fxz - 1i.*fyz).*pi.*(vx + 1i.*vy).*vz;
                    case 2
                        f = -(sqrt(2/5).*(fxx - 1i.*fxy - fyy).*pi.*(vx + 1i.*vy).*vz);
                end
            case 2
                switch k
                    case -2
                        f = ((fxx + 1i.*fxy - fyy).*pi.*(vx + 1i.*vy).^2)/sqrt(10);
                    case -1
                        f = ((fxz + 1i.*fyz).*pi.*(vx + 1i.*vy).^2)/sqrt(10);
                    case 0
                        f = -(((fxx + fyy - 2.*fzz).*pi.*(vx + 1i.*vy).^2)/sqrt(15));
                    case 1
                        f = -(((fxz - 1i.*fyz).*pi.*(vx + 1i.*vy).^2)/sqrt(10));
                    case 2
                        f = ((fxx - 1i.*fxy - fyy).*pi.*(vx + 1i.*vy).^2)/sqrt(10);
                end
        end
end

% drot = 1/20/6;
% d0 = 0;
% d1 = (1/20-1/40)/6;
% t = 0:100;
% zpp0 = 0*t;
% zpp1 = zpp0;
% zpo0 = zpp0;
% zpo1 = zpp0;
% zop0 = zpp0;
% zop1 = zpp0;
% for L=0:6
%     for M=-L:L
%         zpp0 = zpp0 + real(pp(L+1,abs(M)+1))*exp(-L*(L+1)*drot*t - M^2*d0*t);
%         zpo0 = zpo0 + real(po(L+1,abs(M)+1))*exp(-L*(L+1)*drot*t - M^2*d0*t);
%         zop0 = zop0 + real(op(L+1,abs(M)+1))*exp(-L*(L+1)*drot*t - M^2*d0*t);        
% 
%         zpp1 = zpp1 + real(pp(L+1,abs(M)+1))*exp(-L*(L+1)*drot*t - M^2*d1*t);
%         zpo1 = zpo1 + real(po(L+1,abs(M)+1))*exp(-L*(L+1)*drot*t - M^2*d1*t);
%         zop1 = zop1 + real(op(L+1,abs(M)+1))*exp(-L*(L+1)*drot*t - M^2*d1*t);        
%     end
% end
% plot(t,zpp0,t,zpp1,t,zpo0,t,zpo1,t,zop0,t,zop1)
% xlabel('time (ns)');
% ylabel('correlation (a.u.)')

