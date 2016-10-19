% Program FCS for ab initio calculation of FCS data
close all

if iso % fast or slow rotational diffusion
    nv = gold(wavelength==635); dv=0.05; nv = [nv 1.5]; 
    [exc, rho, z] = NewEXC([0 0.1], [focpos-1 focpos+5], NA, n0, n1, n1, [], focpos+5, [], lamex, over, focpos, atf);
    delta = (z(1,exc(1,:)==max(exc(1,:)))-focpos)/2
    [exc, rho, z] = NewEXC(rhofield, zfield, NA, n0, n1, nv, [], d, dv, lamex, over, focpos-delta, atf);
    nv = gold(wavelength==670); nv = [nv 1.5]; 
    cef = NewCEF(rho, z, NA, n0, n1, nv, [], d, dv, lamem, mag, focpos-delta, av, zpin, atf, 1);
    f0 = exc.*cef;
else
    [exc, rho, z, fx0 fx2, fz] = NewEXC(rhofield, zfield, NA, n0, n1, nv, [], d, dv, lamex, over, focpos, atf);
    [cef, rp, rv] = NewCEF(rho, z, NA, n0, n1, nv, [], d, dv, lamem, mag, focpos, av, zpin, atf, 0);
    [lvd,lvu,lpd,lpu,qvd,qvu,qpd,qpu] = LifetimeL(2*pi*z(1,:)/lamem, n0, n1, nv, [], 2*pi*d/lamem, 2*pi*dv/lamem);
    qv = qvd+qvu;
    qp = qpd+qpu;
    qd = qp-qv;
    fx = (abs(fx0).^2 + abs(fx2).^2);
    fz = abs(fz).^2;
    clear f0
    for j=1:length(z(1,:)) 
        if abs(qd(j))>1e-5
            f0(:,j) = (fx(:,j).*((2*rp(:,j)+rv(:,j)).*qp(j) - (5*rp(:,j)-2*rv(:,j)).*qv(j)) + ...
                fz(:,j).*(rv(:,j).*(qv(j)-4*qp(j)) + rp(:,j).*(qp(j)+2*qv(j))))./qd(j).^2/3;            
            f0(:,j) = f0(:,j) + (fz(:,j).*qp(j)-fx(:,j).*qv(j)).*(qp(j).*rv(:,j)-qv(j).*rp(:,j)).*atanh(sqrt(qd(j)./qp(j)))./sqrt(qp(j).*qd(j).^5);
        else
            f0(:,j) = (fx(:,j).*(8*rp(:,j)+2*rv(:,j)) + fz(:,j).*(2*rp(:,j)+3*rv(:,j)))/15/qp(j);
        end
    end
end
    
dif = 5e-5; % diffusion constant in mum^2/time unit
rhom = rho;

Nsub = 3;
lam = 2^(1/Nsub);
timewin = 1e6;
Ncasc = ceil(log2(timewin));
autotime = lam.^(0:Ncasc*Nsub-1)';
modres = [];
dz = z(1,2)-z(1,1); drho = rho(2,1)-rho(1,1); 
f = f0;
z = ones(size(z,2),1)*z(1,:);
rho = rho(:,1)*ones(1,size(rho,1));
% integrating the diffusion equation:
for j = 1:length(autotime)
    tmp = (erf((z-z'+dz/2)/sqrt(4*dif*autotime(j))) - erf((z-z'-dz/2)/sqrt(4*dif*autotime(j))))/2;
    tmp = tmp + (erf((z+2*focpos+z'+dz/2)/sqrt(4*dif*autotime(j))) - erf((z+2*focpos+z'-dz/2)/sqrt(4*dif*autotime(j))))/2;
    tmp = tmp + (erf((z-2*(d-focpos)+z'+dz/2)/sqrt(4*dif*autotime(j))) - erf((z-2*(d-focpos)+z'-dz/2)/sqrt(4*dif*autotime(j))))/2;   
    ss = sum(tmp);
    tmp(:,ss>1) = tmp(:,ss>1)*diag(1./ss(ss>1)); 
    f = f0*tmp;
    
    rr = rho(:,1)*rho(:,1)'/2/dif/autotime(j);
    tmp = (rho.^2+rho'.^2)/4/dif/autotime(j);
    tmp1 = tmp - rr;
    tmp(rr<=700) = besselj(0, i*rr(rr<=700)).*exp(-tmp(rr<=700));
    tmp(rr>700) = exp(-tmp1(rr>700))./sqrt(rr(rr>700)*2*pi);
    tmp = (ones(size(rho(:,1)))*rho(:,1)').*tmp/2/dif/autotime(j)*drho;
    ss = sum(tmp')';
    tmp(ss>1,:) = diag(1./ss(ss>1))*tmp(ss>1,:); 
    f = tmp*f;
    modres = [modres; sum(sum(rhom.*f.*f0))];

    %plot(z(1,:),rho(:,1)'*f); drawnow
    %pcolor(rho(:,1)*ones(1,size(z,2)),ones(size(rho,1),1)*z(1,:),f); shading interp; axis image; drawnow
    semilogx(autotime(1:j), modres/modres(1)); axis([autotime([1 end])' 0 1]); drawnow;
end
z = ones(size(rho,1),1)*z(1,:);
rho = rho(:,1)*ones(1,size(z,2));
modres = modres/modres(1);
dz = z(1,2)-z(1,1); drho = rho(2,1)-rho(1,1); 
veff = 2*pi*drho*dz*sum(sum(rho.*f0))^2/sum(sum(rho.*f0.^2));
semilogx(autotime,modres);
