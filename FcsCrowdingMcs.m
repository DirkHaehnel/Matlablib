% Program for modeling FCS in crowded environment

close all

boxlen = 1;

if ~exist('u','var')
    lamex = 0.532; % excitation wavelength
    resolution = [lamex/0.01 lamex/0.05]; % grid resoltuion for MDF calculation
    rhomax = 3*boxlen; % max radius of considered volume
    zmax = 5*boxlen; % max z-extension of considered volume
    rhofield =  [-(lamex/resolution(1))/2 sqrt(2)*rhomax + (lamex/resolution(1))/2];
    focpos = 45; % focus position in solution (arbitrary number)
    zfield = focpos + [-(lamex/resolution(2))/2 zmax + (lamex/resolution(2))/2];
    NA = 1.2; % numerical aperture
    n0 = 1.33; % ref index of solution/immersion medium
    n = n0;
    n1 = n0;
    d0 = [];
    d = 0;
    d1 = [];
    mag = 60; % magnification
    fd = 164.5e3/mag; % focal distance of objective
    over = 5e3; % radius of laser beam
    ring = [];
    av = 70/2; % radius of confocal aperture in mum
    lamem = 0.570; % center emission wavelength
    zpin = 0; % pinhole position

    exc = GaussExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, [], resolution);
    mdf = GaussExc2MDF(exc, NA, n0, n, n1, focpos, lamem, mag, av, zpin);
    [modres, autotime] = FCS([exc.rho(:,1:end-1), exc.rho],[-exc.z(:,end:-1:2)+focpos exc.z-focpos],[mdf.volx(:,end:-1:2,:), mdf.volx],[mdf.voly(:,end:-1:2,:), mdf.voly]);

    u = squeeze(mdf.volx(:,:,1) + mdf.voly(:,:,1));
    close
end

% physics
BoltzmannConstant = 1.380650324e-23*1e7;
AvogadroConstant = 6.0221419947e23;
Temperature = 273.15 + 36;
BeadDiam = 0.15; % bead diam in mum
Sigma = 0.05; % bead boundary broadening in mum
Eta = 1e-2; % solvent viscosity dyn/cm^2*s

% numerics
dt = 1e-5; % time step in seconds
FillFactor = 0.5; % filling factor of large particles
K = round(FillFactor*boxlen^3/(4/3*pi*(BeadDiam/2)^3)*(pi/sqrt(18))); % number of molecules in box
Kt = 100; % number of tracers
N = 2^17; % max number of time channels
time = (1:N)*dt;
dif = BoltzmannConstant*Temperature/(3*pi*Eta*BeadDiam*1e-4)*1e4^2; % diffusion constant in mum^2/time unit
step = sqrt(2*dif*dt);

% generate force look-up table
dr = 1e-3*boxlen;
r = 0:dr:(2*boxlen*sqrt(3)+dr);
f = (BeadDiam-r)/2; f(f<0)=0; % hard core
f = conv([f(end:-1:2) f],exp(-[-r(end:-1:2) r].^2/2/Sigma^2)/sum(exp(-[-r(end:-1:2) r].^2/2/Sigma^2))); % Gaussian smoothing
f = f((end+1)/2:end-length(r)+1);
f = -f./r;
f(1) = 0;
ft = 0.2*f; % force on tracer (e.g. 20 % of inter-bead force)

% initialization
x = 2*rhomax*(rand(1,K)-0.5);
y = 2*rhomax*(rand(1,K)-0.5);
z = 2*zmax*(rand(1,K)-0.5);
xt = 2*rhomax*(rand(1,Kt)-0.5);
yt = 2*rhomax*(rand(1,Kt)-0.5);
zt = 2*zmax*(rand(1,Kt)-0.5);
fluo = zeros(N,Kt);
col = ones(27*K,1); row = ones(1,K); rowt = ones(1,Kt);
rh = r(r<=boxlen)';
hh = 0*rh;

for jt = 1:N
    tic
     
    % interaction
    xx = mod(x',boxlen);
    yy = mod(y',boxlen);
    zz = mod(z',boxlen);
    xx = repmat([xx-boxlen; xx; xx+boxlen],9,1);
    yy = repmat([repmat(yy-boxlen,3,1); repmat(yy,3,1); repmat(yy+boxlen,3,1)],3,1);
    zz = [repmat(zz-boxlen,9,1); repmat(zz,9,1); repmat(zz+boxlen,9,1)];
    xxt = xx*rowt-col*mod(xt,boxlen);
    xx = xx*row-col*mod(x,boxlen);
    yyt = yy*rowt-col*mod(yt,boxlen);
    yy = yy*row-col*mod(y,boxlen);
    zzt = zz*rowt-col*mod(zt,boxlen);
    zz = zz*row-col*mod(z,boxlen);

    % motion
    rr = sqrt(xx.^2 + yy.^2 + zz.^2);
    hh = hh + mHist(rr(:),rh); % pair correlation function
    % plot(rh(2:end-1),hh(2:end-1)./rh(2:end-1).^2/(4*pi*dr)/K/jt,rh(2:end-1),conc*ones(size(rh(2:end-1))),':'); drawnow
    rr = f(1+floor(rr/dr));
    x = mod(rhomax + x + step*randn(1,K) + sum(xx.*rr), 2*rhomax) - rhomax;
    y = mod(rhomax + y + step*randn(1,K) + sum(yy.*rr), 2*rhomax) - rhomax;
    z = mod(zmax + z + step*randn(1,K) + sum(zz.*rr), 2*zmax) - zmax;
    rr = sqrt(xxt.^2 + yyt.^2 + zzt.^2);
    rr = ft(1+floor(rr/dr));
    xt = mod(rhomax + xt + step*randn(1,Kt) + sum(xxt.*rr), 2*rhomax) - rhomax;
    yt = mod(rhomax + yt + step*randn(1,Kt) + sum(yyt.*rr), 2*rhomax) - rhomax;
    zt = mod(zmax + zt + step*randn(1,Kt) + sum(zzt.*rr), 2*zmax) - zmax;

    % fluorescence
    rt = sqrt(xt.^2 + yt.^2);
    fluo(jt,:) = interp2(exc.z-focpos, exc.rho, u, abs(zt), rt);
    fluo(jt,rt>exc.rho(end,1) | abs(zt)>exc.z(1,end)) = 0;

    toc
end

flc = fft(fluo);
flc = abs(ifft(flc.*conj(flc)));
auto = sum(flc(1:N,:),2);
time = time(2:end/2);
auto = auto(2:end/2);
semilogx(time,auto/max(auto),autotime*5e-5/dif,modres)
xlabel('time (s)'); ylabel('autocorrelation (a.u.)')

figure
rh = rh(2:end-1);
hh = hh(2:end-1)./rh.^2/(4*pi*dr)/K/N;
plot(rh,hh,rh,conc*ones(size(rh)),':');
xlabel('radius (\mum)'); ylabel('pair correlation function ( \mum^{-3})')

