% program for the integration of coupled Nernst-Planck equations 
% in cylindrical coordinates; all units Centimeter-Gram-Seconds

close all
clear all

SoL = 299792458;
Volt = 1e6/SoL;
Coulomb = 10*SoL;
BoltzmannConstant = 1.380650324e-23*1e7;
AvogadroConstant = 6.0221419947e23;
ElectronCharge = 1.60217646263e-19*Coulomb; % in statC
FaradayConstant = AvogadroConstant*ElectronCharge;
Temperature = 273.15 + 36;
DielectricConstant = 80;

Tmax = 0.1; % max integration time in sec
rad = 10e-4; % radius of cell

% diffusion coefficients
DifK = 1.96e-5; 
DifNa = 1.33e-5; 
DifCl = 2.02e-5; 

% bulk concentrations in mol/cm^3
cK0 = 5e-3*1e-3; 
cNa0 = 140e-3*1e-3; 
cCl0 = 145e-3*1e-3; 

% K influx
kappa = 1e-9; % M/s

% initial conditions
cK = 0; 
cNa = 0; 
cCl = 0; 

% time step
dt = 1e-3; % time step in sec

% relaxation matrix
fac = 4*pi*AvogadroConstant*ElectronCharge^2/BoltzmannConstant/Temperature/DielectricConstant;
M = fac*[-DifK*cK0*[1 1 -1]; -DifNa*cNa0*[1 1 -1]; DifCl*cCl0*[1 1 -1]];
M = expm(dt*M)';

% rescaling 
DifK = 2*DifK*dt/rad^2; 
DifNa = 2*DifNa*dt/rad^2; 
DifCl = 2*DifCl*dt/rad^2; 
kappa = kappa*dt;
t = 0;

% start of integration
resK = 0;
resCl = 0;
resNa = 0;
for j=1:ceil(Tmax/dt)
    current = 0;
    % K-influx:
    cK = cK + kappa;

    % Charge relaxation
    tmp = [cK cNa cCl]*M;
    cK = tmp(1);
    cNa = tmp(2);
    cCl = tmp(3);

    % Diffusion
    cK = cK/(1+DifK);
    cNa = cNa/(1+DifNa);
    cCl = cCl/(1+DifCl);

    t(end+1) = t(end) + dt;
    resK(end+1) = cK;
    resNa(end+1) = cNa;
    resCl(end+1) = cCl;
    plot(t,resK*1e12,t,resNa*1e12,t,resCl*1e12);
    xlabel('time [s]')
    ylabel('concentrations [nM]')
    ax = axis;
    axis([0 Tmax ax(3) ax(4)])
    drawnow
end

