% ParaboloidPerfect Raytracing

close all
clear all

load metals
%nal = silver(wavelength==670);
nal = 1;

nm = 1.33;
ng = 1.51;
zv = pi/1000;
rv = 0*pi;

a = 1/(1e4*2*pi); % paraboloid parameter
zf = -3e5; % lens position
lfoc = 1e6; % lens focal length

theta_c = asin(nm/ng);
theta = (theta_c:pi/1e6:pi/2)'; % polar angle of emission; only SAF is considered
%theta = (pi/6:pi/2e3:theta_c)'; % polar angle of emission; only subcritical emission is considered
dt = mean(diff(theta));

phiv = (0:pi/180:pi)'; % azimuthal angle of emission
zz = [zeros(length(theta),2) ones(length(theta),1)]; % unit vector along z

[x,y] = meshgrid(-200:200,0:200); % detector plane

nmax = 200; % sets max Fourier components
nv = (-nmax:nmax)/nmax;

[v,pc,ps] = DipoleL(theta,zv,ng,nm,nm,[],zv,[]);

phase = unwrap(angle(v));
grad = gradient(phase,dt);
z0 = -grad./sin(theta)/ng;
z0(1) = z0(2);
z0(end) = z0(end-1);

w = sin(theta/2);
w = sqrt(nal^2-ng^2*(1-w.^2))./w/ng;
s_phase = -2*atan(imag(w));
p_phase = -2*atan(ng^2*imag(w));

rho0 = z0.*tan(theta);
rho0prime = gradient(rho0,dt);
p_phaseprime = gradient(p_phase,dt);

rho = 1/a*ones(length(theta),1);
for j=length(theta):-1:2
    rho(j-1) = rho(j) + dt*(rho0prime(j)*cos(theta(j))+(rho(j)-rho0(j))/sin(theta(j)));
    %rho(j-1) = rho(j) + dt*(rho(j)/sin(theta(j)));
end
z = -(rho-rho0)./tan(theta);
plot(rho,z,rho,1/2/a*(1-a^2*rho.^2),[0 max(rho)],[0 -max(rho)/tan(theta_c)],':')


% in plane dipole 
fxx = zeros(2*nmax+1,2*nmax+1); fyx = fxx; gxx = fxx; gyx = fxx;
fxy = fxx; fyy = fxx; gxy = fxx; gyy = fxx;
phase = unwrap(angle(ps));
grad = gradient(phase)/mean(diff(theta));
z0 = -grad./sin(theta)/ng;
z0(1) = z0(2);
z0(end) = z0(end-1);
z0(z0>100) = min(z0);


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
    
    px = tp.*repmat(pc*cos(phi),1,3) + ts.*repmat(ps*sin(phi),1,3);
    py = tp.*repmat(pc*sin(phi),1,3) - ts.*repmat(ps*cos(phi),1,3);

    % reflection by the paraboloid:
    
    lam = -(a*sin(theta).*(rr(:,1)*cos(phi)+rr(:,2)*sin(phi))-cos(theta))./(a*sin(theta).^2);
    lam = lam + sqrt(lam.^2 + (1-a^2*rr(:,1).^2)./(a.^2*sin(theta).^2));
    
    tst = ~(imag(lam)==0);
    rr(tst,:) = [];
    nn(tst,:) = [];
    px(tst,:) = [];
    py(tst,:) = [];
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
    w = sqrt(nal^2-ng^2*(1-w.^2))./w/ng;
    s_phase = exp(-i*2*atan(imag(w)));
    p_phase = exp(-i*2*atan(ng^2*imag(w)));
    nn = nn - 2*nvec.*(sum(nn.*nvec,2)*[1 1 1]);
    tp_prime = -cross(ts,nn);
    px = ts.*((sum(ts.*px,2).*s_phase)*[1 1 1]) + tp_prime.*((sum(tp.*px,2).*p_phase)*[1 1 1]);
    py = ts.*((sum(ts.*py,2).*s_phase)*[1 1 1]) + tp_prime.*((sum(tp.*py,2).*p_phase)*[1 1 1]);
    
    % focusing by ideal lens
    
    lam = (zf-rr(:,3))./nn(:,3);
    rr = rr + repmat(lam,1,3).*nn;
    
	%plot3(rr(:,1),rr(:,2),rr(:,3));
    
    ts = cross(nn,zz);
    ts = ts./(sqrt(sum(ts.*ts,2))*[1 1 1]);
    tp = cross(ts,nn);    
    rho = sqrt(sum(rr(:,1:2).^2,2));
    nrho = sum(rr(:,1:2).*nn(:,1:2),2)./sqrt(sum(rr(:,1:2).^2,2));
    nphi = real(sqrt(sum(nn(:,1:2).^2,2)-nrho.^2));
    sprime = (1/lfoc - nrho./rho.*nn(:,3));
    bend = atan(rho.*sprime);
    tmp = sqrt(1 - nphi.^2);
    nrho = -tmp.*sin(bend); nz = -tmp.*cos(bend);
    nn = [(rr(:,1).*nrho - rr(:,2).*nphi)./sqrt(sum(rr(:,1:2).^2,2)) ...
        (rr(:,2).*nrho + rr(:,1).*nphi)./sqrt(sum(rr(:,1:2).^2,2)) nz];
    tp_prime = cross(ts,nn);
    px = ts.*(sum(ts.*px,2)*[1 1 1]) + tp_prime.*(sum(tp.*px,2)*[1 1 1]);
    py = ts.*(sum(ts.*py,2)*[1 1 1]) + tp_prime.*(sum(tp.*py,2)*[1 1 1]);
    
    % detector plane
    
    lam = -lfoc./nn(:,3);
    rr = rr + repmat(lam,1,3).*nn;
    
	%plot3(rr(:,1),rr(:,2),rr(:,3));
    %hold off
    
    tmp = exp(-i*(nn(:,1).*rr(:,1) + nn(:,2).*rr(:,2)));
    fxx = fxx + mhist2(nn(:,1),nn(:,2),nv,nv,px(:,1).*tmp);
    fyx = fyx + mhist2(nn(:,1),nn(:,2),nv,nv,px(:,2).*tmp);
    fxy = fxy + mhist2(nn(:,1),nn(:,2),nv,nv,py(:,1).*tmp);
    fyy = fyy + mhist2(nn(:,1),nn(:,2),nv,nv,py(:,2).*tmp);
	if  ~(phi==0 || phi==pi)
        fxx = fxx + mhist2(nn(:,1),-nn(:,2),nv,nv,px(:,1).*tmp);
        fyx = fyx + mhist2(nn(:,1),-nn(:,2),nv,nv,-px(:,2).*tmp);
        fxy = fxy + mhist2(nn(:,1),-nn(:,2),nv,nv,-py(:,1).*tmp);
        fyy = fyy + mhist2(nn(:,1),-nn(:,2),nv,nv,py(:,2).*tmp);
    end
    px = cross(nn,px);
    py = cross(nn,py);
    gxx = gxx + mhist2(nn(:,1),nn(:,2),nv,nv,px(:,1).*tmp);
    gyx = gyx + mhist2(nn(:,1),nn(:,2),nv,nv,px(:,2).*tmp);
    gxy = gxy + mhist2(nn(:,1),nn(:,2),nv,nv,py(:,1).*tmp);
    gyy = gyy + mhist2(nn(:,1),nn(:,2),nv,nv,py(:,2).*tmp);
	if  ~(phi==0 || phi==pi)    
        gxx = gxx + mhist2(nn(:,1),-nn(:,2),nv,nv,-px(:,1).*tmp);
        gyx = gyx + mhist2(nn(:,1),-nn(:,2),nv,nv,px(:,2).*tmp);
        gxy = gxy + mhist2(nn(:,1),-nn(:,2),nv,nv,py(:,1).*tmp);
        gyy = gyy + mhist2(nn(:,1),-nn(:,2),nv,nv,-py(:,2).*tmp);
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

% out of plane dipole 
fxz = zeros(2*nmax+1,2*nmax+1); fyz = fxz; gxz = fxz; gyz = fxz;
phase = unwrap(angle(v));
grad = gradient(phase)/mean(diff(theta));
z0 = -grad./sin(theta)/ng;
z0(1) = z0(2);
z0(end) = z0(end-1);
z0(z0>100) = min(z0);

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
    
    pz = tp.*repmat(v,1,3);

    % reflection by the paraboloid:
    
    lam = -(a*sin(theta).*(rr(:,1)*cos(phi)+rr(:,2)*sin(phi))-cos(theta))./(a*sin(theta).^2);
    lam = lam + sqrt(lam.^2 + (1-a^2*rr(:,1).^2)./(a.^2*sin(theta).^2));
    
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
    w = sqrt(nal^2-ng^2*(1-w.^2))./w/ng;
    s_phase = exp(-i*2*atan(imag(w)));
    p_phase = exp(-i*2*atan(ng^2*imag(w)));
    nn = nn - 2*nvec.*(sum(nn.*nvec,2)*[1 1 1]);
    tp_prime = -cross(ts,nn);
    pz = ts.*((sum(ts.*pz,2).*s_phase)*[1 1 1]) + tp_prime.*((sum(tp.*pz,2).*p_phase)*[1 1 1]);
    
    % focusing by ideal lens
    
    lam = (zf-rr(:,3))./nn(:,3);
    rr = rr + repmat(lam,1,3).*nn;
    
	%plot3(rr(:,1),rr(:,2),rr(:,3));
    
    ts = cross(nn,zz);
    ts = ts./(sqrt(sum(ts.*ts,2))*[1 1 1]);
    tp = cross(ts,nn);    
    rho = sqrt(sum(rr(:,1:2).^2,2));
    nrho = sum(rr(:,1:2).*nn(:,1:2),2)./sqrt(sum(rr(:,1:2).^2,2));
    nphi = real(sqrt(sum(nn(:,1:2).^2,2)-nrho.^2));
    sprime = (1/lfoc - nrho./rho.*nn(:,3));
    bend = atan(rho.*sprime);
    tmp = sqrt(1 - nphi.^2);
    nrho = -tmp.*sin(bend); nz = -tmp.*cos(bend);
    nn = [(rr(:,1).*nrho - rr(:,2).*nphi)./sqrt(sum(rr(:,1:2).^2,2)) ...
        (rr(:,2).*nrho + rr(:,1).*nphi)./sqrt(sum(rr(:,1:2).^2,2)) nz];
    tp_prime = cross(ts,nn);
    pz = ts.*(sum(ts.*pz,2)*[1 1 1]) + tp_prime.*(sum(tp.*pz,2)*[1 1 1]);
    
    % detector plane
    
    lam = -lfoc./nn(:,3);
    rr = rr + repmat(lam,1,3).*nn;
    
	%plot3(rr(:,1),rr(:,2),rr(:,3));
    %hold off
    
    tmp = exp(-i*(nn(:,1).*rr(:,1) + nn(:,2).*rr(:,2)));
    fxz = fxz + mhist2(nn(:,1),nn(:,2),nv,nv,pz(:,1).*tmp);
    fyz = fyz + mhist2(nn(:,1),nn(:,2),nv,nv,pz(:,2).*tmp);
	if  ~(phi==0 || phi==pi)
        fxz = fxz + mhist2(nn(:,1),-nn(:,2),nv,nv,pz(:,1).*tmp);
        fyz = fyz + mhist2(nn(:,1),-nn(:,2),nv,nv,-pz(:,2).*tmp);
    end
    pz = cross(nn,pz);
    gxz = gxz + mhist2(nn(:,1),nn(:,2),nv,nv,pz(:,1).*tmp);
    gyz = gyz + mhist2(nn(:,1),nn(:,2),nv,nv,pz(:,2).*tmp);
	if  ~(phi==0 || phi==pi)    
        gxz = gxz + mhist2(nn(:,1),-nn(:,2),nv,nv,-pz(:,1).*tmp);
        gyz = gyz + mhist2(nn(:,1),-nn(:,2),nv,nv,pz(:,2).*tmp);
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
