% ParaboloidMDF Raytracing

close all
%clear all
raytracetrack = false;

% load metals
% nal = gold(wavelength==670);
nal = 1;

nm = 1; 
ng = 1.52;
zv = pi/1000;
rv = 0*pi;

a = 1/(1e4*2*pi); % paraboloid parameter
h = 0*pi; % paraboloid defocusing
zf = -5e5; % lens position
lfoc = 1e6; % lens focal length
%focpos = 0e2; % detector defocusing

theta_c = asin(nm/ng); % TIRF angle
%theta = (theta_c:pi/4e3:pi/2)'; % polar angle of emission; only SAF is considered
theta = (pi/6:pi/4e3:theta_c)'; % polar angle of emission; only subcritical emission is considered
%theta = (pi/6:pi/4e3:pi/2)'; % polar angle of emission; both subcritical and SAF 

phiv = (0:pi/180:pi)'; % azimuthal angle of emission
zz = [zeros(length(theta),2) ones(length(theta),1)]; % unit vector along z

tmp = 150;
%tmp = 200;
[x,y] = meshgrid((-tmp:tmp)+rv*lfoc/zf,0:tmp); % detector plane

nmax = 200; % sets max Fourier components
nv = (-nmax:nmax)/nmax;

[v,pc,ps] = DipoleL(theta,zv,ng,nm,nm,[],zv,[]);

% in plane dipole 
fxx = zeros(2*nmax+1,2*nmax+1); fyx = fxx; gxx = fxx; gyx = fxx;
fxy = fxx; fyy = fxx; gxy = fxx; gyy = fxx;
phase = unwrap(angle(ps));
grad = gradient(phase)/mean(diff(theta));
z0 = -grad./sin(theta)/ng;
z0(1) = z0(2);
z0(end) = z0(end-1);

tmp = exp(1i*ng*z0.*sin(theta).^2./cos(theta));
pc = pc.*tmp;
ps = ps.*tmp;

% semilogy(theta, z0.*tan(theta)/2/pi*0.67)
% semilogy(theta, z0/2/pi*0.67)

tic
for j=1:length(phiv)
    phi = phiv(j);
    
    rr = [rv*ones(length(theta),1)+z0.*tan(theta)*cos(phi) z0.*tan(theta)*sin(phi) zeros(length(theta),1)];
    nn = [sin(theta)*cos(phi) sin(theta)*sin(phi) -cos(theta)];
    
    if raytracetrack
        pos1 = rr;
    end
    
    ts = cross(nn,zz);
    ts = ts./(sqrt(sum(ts.*ts,2))*[1 1 1]);
    tp = cross(ts,nn);

    tmp = sin(theta); % weighting with sin(theta)
    px = tp.*repmat(tmp.*pc*cos(phi),1,3) + ts.*repmat(tmp.*ps*sin(phi),1,3);
    py = tp.*repmat(tmp.*pc*sin(phi),1,3) - ts.*repmat(tmp.*ps*cos(phi),1,3);

    % reflection by the paraboloid:
    
    lam = -(a*sin(theta).*(rr(:,1)*cos(phi)+rr(:,2)*sin(phi))-cos(theta))./(a*sin(theta).^2);
    lam = lam + sqrt(lam.^2 + (1+2*a*h-a^2*rr(:,1).^2)./(a.^2*sin(theta).^2));
    
    tst = ~(imag(lam)==0);
    rr(tst,:) = [];
    nn(tst,:) = [];
    px(tst,:) = [];
    py(tst,:) = [];
    lam(tst) = [];
    
    pos1(tst) = 0;
    
    rr = rr + repmat(lam,1,3).*nn;

    if raytracetrack
        pos2 = rr;
    end
    
    %plot3(rr(:,1),rr(:,2),rr(:,3));
    %hold on
    
    %plot(rr(:,1),rr(:,3),'-o',1/a:100:2e5,1/2/a*(1-a^2*(1/a:100:2e5).^2))
    
    nvec = [a*rr(:,1), a*rr(:,2), ones(size(rr,1),1)]./(sqrt(1+a^2*sum(rr(:,1:2).^2,2))*[1 1 1]);
    ts = cross(nn,nvec);
    ts = ts./(sqrt(sum(ts.*ts,2))*[1 1 1]);
    tp = cross(ts,nn);
    w = sum(nn.*nvec,2);
    [rp, rs] = Fresnel(ng*w,ng,nal);
    nn = nn - 2*nvec.*(sum(nn.*nvec,2)*[1 1 1]);
    tp_prime = cross(ts,nn);
    px = ts.*((sum(ts.*px,2).*rs.')*[1 1 1]) + tp_prime.*((sum(tp.*px,2).*rp.')*[1 1 1]);
    py = ts.*((sum(ts.*py,2).*rs.')*[1 1 1]) + tp_prime.*((sum(tp.*py,2).*rp.')*[1 1 1]);
    
    % focusing by ideal lens

    lam = (zf-rr(:,3))./nn(:,3);
    rr = rr + repmat(lam,1,3).*nn;

    if raytracetrack
        pos3 = rr;
    end

    %plot3(rr(:,1),rr(:,2),rr(:,3));

    ts = cross(nn,zz);
    tst = sqrt(sum(ts.*ts,2));
    ts(tst<1e-10,1) = 0; ts(tst<1e-10,2) = -1; ts(tst<1e-10,3) = 0;
    ts = ts./(sqrt(sum(ts.*ts,2))*[1 1 1]);
    tp = cross(ts,nn);
    rho = sqrt(sum(rr(:,1:2).^2,2));
    nrho = sum(rr(:,1:2).*nn(:,1:2),2)./rho;
    nphi = real(sqrt(sum(nn(:,1:2).^2,2)-nrho.^2));
    sprime = (1/lfoc - nrho./rho.*nn(:,3));
    bend = atan(rho.*sprime);
    tmp = sqrt(1 - nphi.^2);
    nrho = -tmp.*sin(bend); nz = -tmp.*cos(bend);
    nn = [(rr(:,1).*nrho - rr(:,2).*nphi)./rho (rr(:,2).*nrho + rr(:,1).*nphi)./rho nz];
    tp_prime = cross(ts,nn);
    px = ts.*(sum(ts.*px,2)*[1 1 1]) + tp_prime.*(sum(tp.*px,2)*[1 1 1]);
    py = ts.*(sum(ts.*py,2)*[1 1 1]) + tp_prime.*(sum(tp.*py,2)*[1 1 1]);

    % detector plane

    lam = -(lfoc+focpos)./nn(:,3);
    rr = rr + repmat(lam,1,3).*nn;
    
    %     % focusing by aplanatic objective
    %
    %     lam0 = lam;
    %     lam = (rr(:,3)-focpos)./nn(:,3);
    %     tmp = (exp(i*(ng*lam0-lam))*[1 1 1]);
    %     px = px.*tmp;
    %     py = py.*tmp;
    %     ts = [-nn(:,2), nn(:,1), zeros(size(nn,1),1)]; ts = ts./(sqrt(sum(ts.^2,2))*[1 1 1]);
    %     tp = [-nn(:,1).*nn(:,3), -nn(:,2).*nn(:,3), sum(nn(:,1:2).^2,2)]; tp = tp./(sqrt(sum(tp.^2,2))*[1 1 1]);
    %     rho = sqrt(sum(rr(:,1:2).^2,2));
    %     nrho = sum(rr(:,1:2).*nn(:,1:2),2)./(rho+(rho==0));
    %     nphi = sum([-rr(:,2) rr(:,1)].*nn(:,1:2),2)./(rho+(rho==0));
    %     nrho = nrho/mag;
    %     nz = -sqrt(1 - nrho.^2 - nphi.^2);
    %     tmp = [-(rr(:,1).*nrho - rr(:,2).*nphi)./(rho+(rho==0)) -(rr(:,2).*nrho + rr(:,1).*nphi)./(rho+(rho==0)) nz];
    %     rr = -(rr - repmat(lam,1,3).*nn)*mag;
    %     ww = sqrt(tmp(:,3)./nn(:,3)); % Richards-Wolf factor
    %     nn = tmp;
    %     tp_prime = -cross(ts,nn);
    %     px = ts.*(ww.*sum(ts.*px,2)*[1 1 1]) + tp_prime.*(ww.*sum(tp.*px,2)*[1 1 1]);
    %     py = ts.*(ww.*sum(ts.*py,2)*[1 1 1]) + tp_prime.*(ww.*sum(tp.*py,2)*[1 1 1]);
    
    if raytracetrack
        pos4 = rr;
        plot3([pos1(:,1) pos2(:,1) pos3(:,1) pos4(:,1)]',[pos1(:,2) pos2(:,2) pos3(:,2) pos4(:,2)]',[pos1(:,3) pos2(:,3) pos3(:,3) pos4(:,3)]')
        axis image
        drawnow
    end
    
    tmp = exp(-1i*(nn(:,1).*rr(:,1) + nn(:,2).*rr(:,2)));
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

vx = exp(1i*nv'*x(1,:));
vy = exp(1i*y(:,1)*nv);
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

% vertical dipole 
fxz = zeros(2*nmax+1,2*nmax+1); fyz = fxz; gxz = fxz; gyz = fxz;
phase = unwrap(angle(v));
grad = gradient(phase)/mean(diff(theta));
z0 = -grad./sin(theta)/ng;
z0(1) = z0(2);
z0(end) = z0(end-1);

v = v.*exp(1i*ng*z0.*sin(theta).^2./cos(theta));

% semilogy(theta, z0.*tan(theta)/2/pi*0.67)
% semilogy(theta, z0/2/pi*0.67)

tic
for j=1:length(phiv)
    phi = phiv(j);
    
    rr = [rv*ones(length(theta),1)+z0.*tan(theta)*cos(phi) z0.*tan(theta)*sin(phi) zeros(length(theta),1)];
    nn = [sin(theta)*cos(phi) sin(theta)*sin(phi) -cos(theta)];
    
    ts = cross(nn,zz);
    ts = ts./(sqrt(sum(ts.*ts,2))*[1 1 1]);
    tp = cross(ts,nn);
    
    pz = tp.*repmat(sin(theta).*v,1,3); % weighting with sin(theta)

    % reflection by the paraboloid:
    
    lam = -(a*sin(theta).*(rr(:,1)*cos(phi)+rr(:,2)*sin(phi))-cos(theta))./(a*sin(theta).^2);
    lam = lam + sqrt(lam.^2 + (1+2*a*h-a^2*rr(:,1).^2)./(a.^2*sin(theta).^2));
    
    tst = ~(imag(lam)==0);
    rr(tst,:) = [];
    nn(tst,:) = [];
    pz(tst,:) = [];
    lam(tst) = [];
    
    rr = rr + repmat(lam,1,3).*nn;
    
    %plot3(rr(:,1),rr(:,2),rr(:,3));
    %hold on
    
    %plot(rr(:,1),rr(:,3),'-o',1/a:100:2e5,1/2/a*(1-a^2*(1/a:100:2e5).^2))
    
    nvec = [a*rr(:,1), a*rr(:,2), ones(size(rr,1),1)]./(sqrt(1+a^2*sum(rr(:,1:2).^2,2))*[1 1 1]);
    ts = cross(nn,nvec);
    ts = ts./(sqrt(sum(ts.*ts,2))*[1 1 1]);
    tp = cross(ts,nn);
    w = sum(nn.*nvec,2);
    [rp, rs] = Fresnel(ng*w,ng,nal);
    nn = nn - 2*nvec.*(sum(nn.*nvec,2)*[1 1 1]);
    tp_prime = cross(ts,nn);
    pz = ts.*((sum(ts.*pz,2).*rs.')*[1 1 1]) + tp_prime.*((sum(tp.*pz,2).*rp.')*[1 1 1]);
    
    % focusing by ideal lens

    lam = (zf-rr(:,3))./nn(:,3);
    rr = rr + repmat(lam,1,3).*nn;

    %plot3(rr(:,1),rr(:,2),rr(:,3));

    ts = cross(nn,zz);
    tst = sqrt(sum(ts.*ts,2));
    ts(tst<1e-10,1) = 0; ts(tst<1e-10,2) = -1; ts(tst<1e-10,3) = 0;
    ts = ts./(sqrt(sum(ts.*ts,2))*[1 1 1]);
    tp = cross(ts,nn);
    rho = sqrt(sum(rr(:,1:2).^2,2));
    nrho = sum(rr(:,1:2).*nn(:,1:2),2)./rho;
    nphi = real(sqrt(sum(nn(:,1:2).^2,2)-nrho.^2));
    sprime = (1/lfoc - nrho./rho.*nn(:,3));
    bend = atan(rho.*sprime);
    tmp = sqrt(1 - nphi.^2);
    nrho = -tmp.*sin(bend); nz = -tmp.*cos(bend);
    nn = [(rr(:,1).*nrho - rr(:,2).*nphi)./rho (rr(:,2).*nrho + rr(:,1).*nphi)./rho nz];
    tp_prime = cross(ts,nn);
    pz = ts.*(sum(ts.*pz,2)*[1 1 1]) + tp_prime.*(sum(tp.*pz,2)*[1 1 1]);

    % detector plane

    lam = -(lfoc+focpos)./nn(:,3);
    rr = rr + repmat(lam,1,3).*nn;
    
    %     % focusing by aplanatic objective
    %
    %     lam0 = lam;
    %     lam = (rr(:,3)-focpos)./nn(:,3);
    %     pz = pz.*(exp(i*(ng*lam0-lam))*[1 1 1]);
    %     ts = [-nn(:,2), nn(:,1), zeros(size(nn,1),1)]; ts = ts./(sqrt(sum(ts.^2,2))*[1 1 1]);
    %     tp = [-nn(:,1).*nn(:,3), -nn(:,2).*nn(:,3), sum(nn(:,1:2).^2,2)]; tp = tp./(sqrt(sum(tp.^2,2))*[1 1 1]);
    %     rho = sqrt(sum(rr(:,1:2).^2,2));
    %     nrho = sum(rr(:,1:2).*nn(:,1:2),2)./(rho+(rho==0));
    %     nphi = sum([-rr(:,2) rr(:,1)].*nn(:,1:2),2)./(rho+(rho==0));
    %     nrho = nrho/mag;
    %     nz = -sqrt(1 - nrho.^2 - nphi.^2);
    %     tmp = [-(rr(:,1).*nrho - rr(:,2).*nphi)./(rho+(rho==0)) -(rr(:,2).*nrho + rr(:,1).*nphi)./(rho+(rho==0)) nz];
    %     rr = -(rr - repmat(lam,1,3).*nn)*mag;
    %     ww = sqrt(tmp(:,3)./nn(:,3)); % Richards-Wolf factor
    %     nn = tmp;
    %     tp_prime = -cross(ts,nn);
    %     pz = ts.*(ww.*sum(ts.*pz,2)*[1 1 1]) + tp_prime.*(ww.*sum(tp.*pz,2)*[1 1 1]);
    
	%plot3(rr(:,1),rr(:,2),rr(:,3));
    %hold off
    
    tmp = exp(-1i*(nn(:,1).*rr(:,1) + nn(:,2).*rr(:,2)));
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

% ParaboloidMDF0.mat: zentriertes Molekül, Glasparaboloid, SAF
% ParaboloidMDF0sub.mat: zentriertes Molekül, Glasparaboloid, subkritische Detektion

return

% field distribution along optical axis
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
