% Program for modeling FCS on interacting particles 

close all
boxlen = 5; % size of one MCS box in mum

if ~exist('u','var')
    lamex = 0.543; % excitation wavelength
    resolution = [lamex/0.025 lamex/0.05]; % spatioal grid resolution for MDF calculation 
    rhomax = boxlen; % radial extent of calculation volume in mum
    zmax = 2*boxlen; % axial extent of calculation volume in mum
    rhofield =  [-(lamex/resolution(1))/2 sqrt(2)*rhomax + (lamex/resolution(1))/2];
    focpos = 45; % position of focus in mum
    zfield = focpos + [-(lamex/resolution(2))/2 zmax+(lamex/resolution(2))/2];
    NA = 0.9; % numerical aperture
    n0 = 1.497; % ref. index of solution/immersion medium
    n = n0;
    n1 = n0;
    d0 = [];
    d = 0;
    d1 = [];
    mag = 40; % magnification
    fd = 164.5e3/mag; % focal distance of objective in mum
    over = 5e3; % radius of laser beam in mum
    ring = [];
    av = 70/2; % radius of confocal aperture in mum
    lamem = 0.570; % center emission wavelength
    zpin = 0; % position of confocal aperture

    exc = GaussExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, [], resolution);
    mdf = GaussExc2MDF(exc, NA, n0, n, n1, focpos, lamem, mag, av, zpin);
    [modres, autotime] = FCS([exc.rho(:,1:end-1), exc.rho],[-exc.z(:,end:-1:2)+focpos exc.z-focpos],[mdf.volx(:,end:-1:2,:), mdf.volx],[mdf.voly(:,end:-1:2,:), mdf.voly]);

    u = squeeze(mdf.volx(:,:,1) + mdf.voly(:,:,1));
    close
    
	veff = 2*DetectionVolume(exc.rho,exc.z,mdf.volx,mdf.voly);
end

% physics
BoltzmannConstant = 1.380650324e-23*1e7;
AvogadroConstant = 6.0221419947e23;
ElectronCharge = 1.60217646263e-19*2997924580;
Temperature = 273.15 + 36;
DielectricConstant = 80;
BeadDiam = 5e-2; % bead diam in mum
BeadCharge = 1e3; % bead charge in elementary charge units
Eta = 1e-2; % solvent viscosity dyn/cm^2*s

% numerics
dt = 1e-5; % time step in seconds
force = dt*BeadCharge^2*ElectronCharge^2/(3*pi*Eta*BeadDiam*1e-4)*1e4^2;
conc = 1; % number of particles per mum^3
% IonStrength = 0.5*beadcharge*conc/veff; % molecules/mum^3
% lambdaDH = 1/sqrt(8*pi*ElectronCharge^2/DielectricConstant/BoltzmannConstant/Temperature*IonStrength); % Debye-Hueckel screening length in mum
lambdaDH = 1; % Debye-Hueckel screening length in mum
K = round(boxlen^3*conc); % number of molecules in box
N = 2^17; % max number of time channels
time = (1:N)*dt;
dif = BoltzmannConstant*Temperature/(3*pi*Eta*BeadDiam*1e-4)*1e4^2; % diffusion constant in mum^2/time unit
auto = zeros(N,1);
step = sqrt(2*dif*dt);

% generate force look-up table
dr = 1e-3*boxlen;
r = 0:dr:(2*boxlen*sqrt(3)+dr);
f = force*exp(-r/lambdaDH).*(1./r+1/lambdaDH)./r.^2; % Debye-Hueckel interaction
f(r<BeadDiam) = 0;
f(r>boxlen) = 0;

% initialization
x = 2*rhomax*(rand(1,K)-0.5);
y = 2*rhomax*(rand(1,K)-0.5);
z = 2*zmax*(rand(1,K)-0.5);
fluo = zeros(N,K);
col = ones(27*K,1); row = ones(1,K);
rh = r(r<=boxlen)';
hh = 0*rh;

for k=1:2
    % random path
    if 1
        for jt = 1:N
            %tic

            % interaction
            xx = mod(x',boxlen);
            xx = repmat([xx-boxlen; xx; xx+boxlen],9,1);
            xx = xx*row-col*mod(x,boxlen);
            yy = mod(y',boxlen);
            yy = repmat([repmat(yy-boxlen,3,1); repmat(yy,3,1); repmat(yy+boxlen,3,1)],3,1);
            yy = yy*row-col*mod(y,boxlen);
            zz = mod(z',boxlen);
            zz = [repmat(zz-boxlen,9,1); repmat(zz,9,1); repmat(zz+boxlen,9,1)];
            zz = zz*row-col*mod(z,boxlen);
            rr = sqrt(xx.^2 + yy.^2 + zz.^2);
            hh = hh + mHist(rr(:),rh); % pair correlation function
            %plot(rh(2:end-1),hh(2:end-1)./rh(2:end-1).^2/(4*pi*dr)/K/jt,rh(2:end-1),conc*ones(size(rh(2:end-1))),':'); drawnow
            rr = f(1+floor(rr/dr));

            % motion
            x = mod(rhomax + x + step*randn(1,K) - sum(xx.*rr), 2*rhomax) - rhomax;
            y = mod(rhomax + y + step*randn(1,K) - sum(yy.*rr), 2*rhomax) - rhomax;
            z = mod(zmax + z + step*randn(1,K) - sum(zz.*rr), 2*zmax) - zmax;

            % fluorescence
            rr = sqrt(x.^2 + y.^2);
            fluo(jt,:) = interp2(exc.z-focpos, exc.rho, u, abs(z), rr);
            fluo(jt,rr>exc.rho(end,1) | abs(z)>exc.z(1,end)) = 0;

            %toc
        end
    else % checking the modeling for the case of no interaction
        x = mod(rhomax + cumsum([x(1,:);step*randn(N-1,K)]), 2*rhomax) - rhomax;
        y = mod(rhomax + cumsum([y(1,:);step*randn(N-1,K)]), 2*rhomax) - rhomax;
        z = mod(zmax + cumsum([z(1,:);step*randn(N-1,K)]), 2*zmax) - zmax;
        
        % fluorescence
        rr = sqrt(x.^2 + y.^2);
        fluo = interp2(exc.z-focpos, exc.rho, u, abs(z), rr);
        fluo(rr>exc.rho(end,1) | abs(z)>exc.z(1,end)) = 0;
    end


    flc = fft(sum(fluo,2));
    flc = abs(ifft(flc.*conj(flc)));
    auto = auto + flc(1:N);
    semilogx(time(2:end/2), auto(2:end/2)); drawnow;
    
    disp(k)
end

time = time(2:end/2);
auto = auto(2:end/2);
semilogx(time,auto/max(auto),autotime*5e-5/dif,(modres+K/8/rhomax^2/zmax*veff)/(1+K/8/rhomax^2/zmax*veff))
xlabel('time (s)'); ylabel('autocorrelation (a.u.)')

figure
rh = rh(2:end-1);
hh = hh(2:end-1)./rh.^2/(4*pi*dr)/K/N/k;
plot(rh,hh,rh,conc*ones(size(rh)),':');
xlabel('radius (\mum)'); ylabel('pair correlation function ( \mum^{-3})')

