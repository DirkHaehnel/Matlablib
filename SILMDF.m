% SIL MDF Raytracing

close all
%clear all

nm = 1.;
ng = 1.8;
zv = pi/1000;
rv = 0*pi;

rad = (4.5e3*2*pi); % radius SIL
h = 0*15*2*pi; % SIL decentring
% focpos = 7.5*pi; % objective position
mag = 40; % objective magnification
NA = 0.7; % objective NA

theta = (pi/8e3:pi/4e3:asin(NA))'; % polar angle of emission

phiv = (0:pi/180:pi)'; % azimuthal angle of emission
zz = [zeros(length(theta),2) ones(length(theta),1)]; % unit vector along z

[x,y] = meshgrid((-1e3:2:1e3)-rv*mag*ng,0:2:1e3); % detector plane

nmax = 2e2; % sets max Fourier components
nv = (-nmax:nmax)/nmax/mag;

[v,pc,ps] = DipoleL(theta,zv,ng,nm,nm,[],zv,[]);

% in plane dipole
fxx = zeros(2*nmax+1,2*nmax+1); fyx = fxx; gxx = fxx; gyx = fxx;
fxy = fxx; fyy = fxx; gxy = fxx; gyy = fxx;
phase = unwrap(angle(ps));
grad = gradient(phase)/mean(diff(theta));
z0 = -grad./sin(theta)/ng;
z0(1) = z0(2);
z0(end) = z0(end-1);

% semilogy(theta, z0.*tan(theta)/2/pi*0.67)
% semilogy(theta, z0/2/pi*0.67)

tic
for j=1:length(phiv)
    phi = phiv(j);

    rr = [rv*ones(length(theta),1)+z0.*tan(theta)*cos(phi) z0.*tan(theta)*sin(phi) zeros(length(theta),1)];
    nn = [sin(theta)*cos(phi) sin(theta)*sin(phi) -cos(theta)];

    ts = repmat([sin(phi) -cos(phi) 0],length(theta),1);
    tp = cross(ts,nn);

    tmp = sin(theta); % weighting with sin(theta)
    px = tp.*repmat(tmp.*pc*cos(phi),1,3) + ts.*repmat(tmp.*ps*sin(phi),1,3);
    py = tp.*repmat(tmp.*pc*sin(phi),1,3) - ts.*repmat(tmp.*ps*cos(phi),1,3);
    
    %     % intersection with SIL:

	[rr, nn, px, py, lam0] = Refraction(rr, nn, px, py, rad, h, 1/ng);
    
    %plot(real(px(:,1)),real(px(:,2))); axis image; drawnow; 
    
    % focusing by aplanatic objective
    
    tst = sin(acos(-nn(:,3)))>NA;
    rr(tst,:) = [];
    nn(tst,:) = [];
    px(tst,:) = [];
    py(tst,:) = [];
    lam0(tst) = [];

    lam = (rr(:,3)-focpos)./nn(:,3);
    tmp = (exp(i*(ng*lam0-lam))*[1 1 1]);
    px = px.*tmp;
    py = py.*tmp;    
    ts = [-nn(:,2), nn(:,1), zeros(size(nn,1),1)]; ts = ts./(sqrt(sum(ts.^2,2))*[1 1 1]);
    tp = [-nn(:,1).*nn(:,3), -nn(:,2).*nn(:,3), sum(nn(:,1:2).^2,2)]; tp = tp./(sqrt(sum(tp.^2,2))*[1 1 1]);
    rho = sqrt(sum(rr(:,1:2).^2,2));
    nrho = sum(rr(:,1:2).*nn(:,1:2),2)./(rho+(rho==0));
    nphi = sum([-rr(:,2) rr(:,1)].*nn(:,1:2),2)./(rho+(rho==0));
    nrho = nrho/mag;
    nz = -sqrt(1 - nrho.^2 - nphi.^2);
    tmp = [-(rr(:,1).*nrho - rr(:,2).*nphi)./(rho+(rho==0)) -(rr(:,2).*nrho + rr(:,1).*nphi)./(rho+(rho==0)) nz];
    rr = -(rr - repmat(lam,1,3).*nn)*mag;
    ww = sqrt(tmp(:,3)./nn(:,3)); % Richards-Wolf factor
    nn = tmp;
    tp_prime = -cross(ts,nn);
    px = ts.*(ww.*sum(ts.*px,2)*[1 1 1]) + tp_prime.*(ww.*sum(tp.*px,2)*[1 1 1]);
    py = ts.*(ww.*sum(ts.*py,2)*[1 1 1]) + tp_prime.*(ww.*sum(tp.*py,2)*[1 1 1]);

    %plot3(rr(:,1),rr(:,2),-2*rad*ones(size(rr,1),1));
    %hold off

    %plot(nn(:,1),nn(:,2)); axis image; drawnow; 
    
    tmp = exp(-i*(nn(:,1).*rr(:,1) + nn(:,2).*rr(:,2)));
    fxx = fxx + mHist2(nn(:,1),nn(:,2),nv,nv,px(:,1).*tmp);
    fyx = fyx + mHist2(nn(:,1),nn(:,2),nv,nv,px(:,2).*tmp);
    fxy = fxy + mHist2(nn(:,1),nn(:,2),nv,nv,py(:,1).*tmp);
    fyy = fyy + mHist2(nn(:,1),nn(:,2),nv,nv,py(:,2).*tmp);
    if  ~(phi==0 || phi==pi)
        fxx = fxx + mHist2(nn(:,1),-nn(:,2),nv,nv,px(:,1).*tmp);
        fyx = fyx + mHist2(nn(:,1),-nn(:,2),nv,nv,-px(:,2).*tmp);
        fxy = fxy + mHist2(nn(:,1),-nn(:,2),nv,nv,-py(:,1).*tmp);
        fyy = fyy + mHist2(nn(:,1),-nn(:,2),nv,nv,py(:,2).*tmp);
    end
    px = cross(nn,px);
    py = cross(nn,py);
    gxx = gxx + mHist2(nn(:,1),nn(:,2),nv,nv,px(:,1).*tmp);
    gyx = gyx + mHist2(nn(:,1),nn(:,2),nv,nv,px(:,2).*tmp);
    gxy = gxy + mHist2(nn(:,1),nn(:,2),nv,nv,py(:,1).*tmp);
    gyy = gyy + mHist2(nn(:,1),nn(:,2),nv,nv,py(:,2).*tmp);
    if  ~(phi==0 || phi==pi)
        gxx = gxx + mHist2(nn(:,1),-nn(:,2),nv,nv,-px(:,1).*tmp);
        gyx = gyx + mHist2(nn(:,1),-nn(:,2),nv,nv,px(:,2).*tmp);
        gxy = gxy + mHist2(nn(:,1),-nn(:,2),nv,nv,py(:,1).*tmp);
        gyy = gyy + mHist2(nn(:,1),-nn(:,2),nv,nv,-py(:,2).*tmp);
    end
end
toc

vx = exp(i*nv'*x(1,:));
vy = exp(i*y(:,1)*nv);
exx = vy*fxx*vx;
exy = vy*fxy*vx;
eyx = vy*fyx*vx;
eyy = vy*fyy*vx;
bxx = vy*gxx*vx;
bxy = vy*gxy*vx;
byx = vy*gyx*vx;
byy = vy*gyy*vx;
        
fx = -real(exx.*conj(byx)-eyx.*conj(bxx)); % negative direction of propagation!
fy = -real(exy.*conj(byy)-eyy.*conj(bxy)); % negative direction of propagation!
subplot(131)
pcolor([flipud(x(2:end,:)); x],[-flipud(y(2:end,:)); y],[flipud(fx(2:end,:)); fx]); shading flat; axis image; set(gca,'fontsize',12)
subplot(132)
pcolor([flipud(x(2:end,:)); x],[-flipud(y(2:end,:)); y],[flipud(fy(2:end,:)); fy]); shading flat; axis image; set(gca,'fontsize',12)
drawnow


% out of plane dipole 
fxz = zeros(2*nmax+1,2*nmax+1); fyz = fxz; gxz = fxz; gyz = fxz;
phase = unwrap(angle(v));
grad = gradient(phase)/mean(diff(theta));
z0 = -grad./sin(theta)/ng;
z0(1) = z0(2);
z0(end) = z0(end-1);

tic
for j=1:length(phiv)
    phi = phiv(j);

    rr = [rv*ones(length(theta),1)+z0.*tan(theta)*cos(phi) z0.*tan(theta)*sin(phi) zeros(length(theta),1)];
    nn = [sin(theta)*cos(phi) sin(theta)*sin(phi) -cos(theta)];

    ts = repmat([sin(phi) -cos(phi) 0],length(theta),1);
    tp = cross(ts,nn);

    pz = tp.*repmat(sin(theta).*v,1,3); % weighting with sin(theta)
    
    % intersection with SIL:

	[rr, nn, pz, tmp, lam0] = Refraction(rr, nn, pz, [], rad, h, 1/ng);
    
    %plot(real(px(:,1)),real(px(:,2))); axis image; drawnow; 
    
    % focusing by aplanatic objective
    
    tst = sin(acos(-nn(:,3)))>NA;
    rr(tst,:) = [];
    nn(tst,:) = [];
    pz(tst,:) = [];
    lam0(tst) = [];

    lam = (rr(:,3)-focpos)./nn(:,3);
    pz = pz.*(exp(i*(ng*lam0-lam))*[1 1 1]);
    ts = [-nn(:,2), nn(:,1), zeros(size(nn,1),1)]; ts = ts./(sqrt(sum(ts.^2,2))*[1 1 1]);
    tp = [-nn(:,1).*nn(:,3), -nn(:,2).*nn(:,3), sum(nn(:,1:2).^2,2)]; tp = tp./(sqrt(sum(tp.^2,2))*[1 1 1]);
    rho = sqrt(sum(rr(:,1:2).^2,2));
    nrho = sum(rr(:,1:2).*nn(:,1:2),2)./(rho+(rho==0));
    nphi = sum([-rr(:,2) rr(:,1)].*nn(:,1:2),2)./(rho+(rho==0));
    nrho = nrho/mag;
    nz = -sqrt(1 - nrho.^2 - nphi.^2);
    tmp = [-(rr(:,1).*nrho - rr(:,2).*nphi)./(rho+(rho==0)) -(rr(:,2).*nrho + rr(:,1).*nphi)./(rho+(rho==0)) nz];
    rr = -(rr - repmat(lam,1,3).*nn)*mag;
    ww = sqrt(tmp(:,3)./nn(:,3)); % Richards-Wolf factor
    nn = tmp;
    tp_prime = -cross(ts,nn);
    pz = ts.*(ww.*sum(ts.*pz,2)*[1 1 1]) + tp_prime.*(ww.*sum(tp.*pz,2)*[1 1 1]);

    tmp = exp(-i*(nn(:,1).*rr(:,1) + nn(:,2).*rr(:,2)));
    fxz = fxz + mHist2(nn(:,1),nn(:,2),nv,nv,pz(:,1).*tmp);
    fyz = fyz + mHist2(nn(:,1),nn(:,2),nv,nv,pz(:,2).*tmp);
    if  ~(phi==0 || phi==pi)
        fxz = fxz + mHist2(nn(:,1),-nn(:,2),nv,nv,pz(:,1).*tmp);
        fyz = fyz + mHist2(nn(:,1),-nn(:,2),nv,nv,-pz(:,2).*tmp);
    end
    pz = cross(nn,pz);
    gxz = gxz + mHist2(nn(:,1),nn(:,2),nv,nv,pz(:,1).*tmp);
    gyz = gyz + mHist2(nn(:,1),nn(:,2),nv,nv,pz(:,2).*tmp);
    if  ~(phi==0 || phi==pi)
        gxz = gxz + mHist2(nn(:,1),-nn(:,2),nv,nv,-pz(:,1).*tmp);
        gyz = gyz + mHist2(nn(:,1),-nn(:,2),nv,nv,pz(:,2).*tmp);
    end
end
toc

exz = vy*fxz*vx;
eyz = vy*fyz*vx;
bxz = vy*gxz*vx;
byz = vy*gyz*vx;

fz = -real(exz.*conj(byz)-eyz.*conj(bxz)); % negative direction of propagation!
subplot(133)
pcolor([flipud(x(2:end,:)); x],[-flipud(y(2:end,:)); y],[flipud(fz(2:end,:)); fz]); shading flat; axis image; set(gca,'fontsize',12)
drawnow

colormap hot

return

for jf=1:30 focpos=(jf-7)*0.5*pi; SILMDF; eval(['print -dpng -r300 tst' mint2str(jf,2)]); end

return

if 0
    close
    for j=1:50
        vx = exp(i*nv'*x(1,:));
        vy = exp(i*y(:,1)*nv);
        focpos = -2*pi*j*2e1;
        vz = exp(i*real(sqrt(1-nv'.^2*ones(1,length(nv))-ones(length(nv),1)*nv.^2))*focpos);
        exx = vy*(fxx.*vz)*vx;
        exy = vy*(fxy.*vz)*vx;
        eyx = vy*(fyx.*vz)*vx;
        eyy = vy*(fyy.*vz)*vx;
        bxx = vy*(gxx.*vz)*vx;
        bxy = vy*(gxy.*vz)*vx;
        byx = vy*(gyx.*vz)*vx;
        byy = vy*(gyy.*vz)*vx;

        fx = -real(exx.*conj(byx)-eyx.*conj(bxx)); % negative direction of propagation!
        fy = -real(exy.*conj(byy)-eyy.*conj(bxy)); % negative direction of propagation!
        subplot(131)
        pcolor([flipud(x(2:end,:)); x],[-flipud(y(2:end,:)); y],[flipud(fx(2:end,:)); fx]); shading flat; axis image; set(gca,'fontsize',12)
        subplot(132)
        pcolor([flipud(x(2:end,:)); x],[-flipud(y(2:end,:)); y],[flipud(fy(2:end,:)); fy]); shading flat; axis image; set(gca,'fontsize',12)
        drawnow
    end
end

