function [dist, flux, rhov, zv, delta, conc] = CalciumDiffusion(L, R, rho0, DiffConstant, Charge, IonStrength, InfluxPoint)

% program for the integration of the diffusion equation within a cylindrical 
% potential in cylindrical coordinates
% input variables:
% L - length of cylinder in nm
% R - radius of cylinder in nm
% rho0 - inner radius in nm
% Charge - line charge of GARP
% IonStrength - ionic strength in M
% InfluxPoint - radial position of influx point

nphi = 0; % number of angular components to be considered
nz = 400; % number of axial components to be considered
L = L*1e-7; % length of cylinder in cm
delta = rho0*1e-7; % spatial step size in cm; 
rho0 = round(rho0*1e-7/delta)*delta;
rhov = (rho0 - delta/2):delta:(R*1e-7);
nrho = length(rhov)-2;
zv = 0:1e-8:L;
kaj = pi/L*((0:nz)+0.5);

BoltzmannConstant = 1.380650324e-23*1e7;
AvogadroConstant = 6.0221419947e23;
ElectronCharge = 1.60217646263e-19*2997924580;
Temperature = 273.15 + 36;
DielectricConstant = 80;
ql = Charge*1e7*2*pi*ElectronCharge/DielectricConstant; % statC/cm
ql = ql*2*ElectronCharge/BoltzmannConstant/Temperature;
IonStrength = IonStrength*0.5*AvogadroConstant*1e-3; % molecules/ml
kappa = sqrt(8*pi*ElectronCharge^2/DielectricConstant/BoltzmannConstant/Temperature*IonStrength);
psi = real(i*ql*besselh(0,2,-i*kappa*rhov));
zeta = psi(1);
lam0 = exp(-psi/2).*sqrt(rhov);
lam = -(lam0(1:end-2) + lam0(3:end))./lam0(2:end-1);
M = spdiags(ones(nrho,1),0,nrho,nrho);
for jrho=1:nrho
    if jrho>1 & jrho<nrho
        M(jrho,jrho) = lam(jrho);
        M(jrho,jrho-1) = 1;
        M(jrho,jrho+1) = 1;            
    elseif jrho==nrho
        M(jrho,jrho) = lam(jrho) + sqrt(rhov(end))/sqrt(rhov(end-1))*exp((psi(end-1)-psi(end))/2);
        M(jrho,jrho-1) = 1;
    elseif jrho==1
        M(jrho,jrho) = lam(jrho) + sqrt(rhov(1))/sqrt(rhov(2))*exp((psi(2)-psi(1))/2);;
        M(jrho,jrho+1) = 1;
    end            
end
lam0([1 end]) = [];
psi([1 end]) = [];
rhov([1 end]) = [];

jSource = sum(rhov<=InfluxPoint*1e-7);
conc = zeros(nrho,nz);
conc(jSource,:) = -delta/DiffConstant/pi/L/sqrt(rhov(jSource))*exp(psi(jSource)/2);

for j=1:nz
    conc(:,j) = (M-spdiags(delta^2*kaj(j)^2*ones(nrho,1),0,nrho,nrho))\conc(:,j);
end
conc = conc./((sqrt(rhov).*exp(psi/2))'*ones(1,nz));

dist = zeros(length(rhov),length(zv));
for j=1:nz
    dist = dist + conc(:,j)*cos(kaj(j)*zv);
end
dist(dist<0) = 0;
dist(:,end)=0;
flux = DiffConstant*dist(:,end-1)/mean(diff(zv));

rhov = 1e7*rhov;
zv = 1e7*zv;

% total flux: 2*pi*rhov*flux*1e-7*delta ~ 1

