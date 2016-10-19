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
cyl_thickness = 5e-4; % thickness of enclosing cylinder
cyl_height = 5e-4; % height of enclosing cylinder
cyl_depth = 5e-4; % depth of supporting cylinder (silicon oxide dielectric)
h = 1e-5; % height of gap = FE cell size for integrating the Poisson equation
delta = 2.5e-5; % FE cell size for integrating the Nernst-Planck equations, should be an odd multiple of h

% diffusion coefficients
DifK = 1.96e-5; 
DifNa = 1.33e-5; 
DifCl = 2.02e-5; 

% bulk concentrations in mol/cm^3
cK0 = 5e-3*1e-3; 
cNa0 = 140e-3*1e-3; 
cCl0 = 145e-3*1e-3; 

% inverse membrane conductivity
kappa = 1/1e8*Coulomb/Volt/FaradayConstant/(4*pi*rad^2); % statC/statV/cm^2/s

% internal potential
V0 = 50e-3*1e6/SoL; % statV

jump = round(cyl_thickness/delta);
[z,rho] = meshgrid(delta/2:delta:cyl_height,delta/2:delta:rad+cyl_thickness);
rho = rho(:);
z = z(:);
ind = rho<rad & z>delta;
rho(ind) = []; z(ind) = [];
nn = length(rho);
ind = (1:nn)';
gap = ind(rho<rad);
bottom = ind(z<delta);
outer = ind(rho>(rad+cyl_thickness-delta));
inner = ind(rho<rad+delta & z>delta);
top = ind(z>cyl_height-delta);
nl = ind-1;
nr = ind+1;
nu = ind+jump;
nd = ind-jump;
nl(1) = 1;
nl(inner) = ind(inner);
nr(outer) = ind(outer);
nd(bottom) = ind(bottom);
nu(top) = ind(top);
nu(gap) = ind(gap);
vol = ind(rho>rad); 
rho2 = sqrt(rho);

% Laplace operator for diffusion integrator
L = spdiags(-4*ones(nn,1),0,nn,nn);
for j=1:nn
    L(j,nr(j)) = L(j,nr(j))+1;
    L(j,nl(j)) = L(j,nl(j))+1;
    L(j,nu(j)) = L(j,nu(j))+1;
    L(j,nd(j)) = L(j,nd(j))+1;
end
L(ind(rho>rad & rho<rad+delta & z<delta),ind(rho>rad & rho<rad+delta & z<delta)) = -2 - h/delta;
L(ind(rho>rad & rho<rad+delta & z<delta),ind(rho>rad-delta & rho<rad  & z<delta)) = h/delta;
rho2L = (L*rho2)./rho2;

% boundary conditions
for j=1:length(outer)
    L(outer(j),outer(j)) = L(outer(j),outer(j))-1;
end
for j=1:length(top)
    L(top(j),top(j)) = L(top(j),top(j))-1;
end

% initial conditions
cK = zeros(nn,1); 
cNa = zeros(nn,1); 
cCl = zeros(nn,1); 
qK = zeros(size(gap));

% time step
ddt = 1e-3; % time step in sec
dt = min(0.2*delta^2./[DifK DifNa DifCl]); % integration time step in sec
dt = ddt/floor(ddt/dt);

% relaxation matrix
fac = 4*pi*AvogadroConstant*ElectronCharge^2/BoltzmannConstant/Temperature/DielectricConstant;
M = fac*[-DifK*cK0*[1 1 -1]; -DifNa*cNa0*[1 1 -1]; DifCl*cCl0*[1 1 -1]];
M = expm(dt*M)';

% rescaling 
DifK = DifK/delta^2*dt; 
DifNa = DifNa/delta^2*dt; 
DifCl = DifCl/delta^2*dt; 
kappa = kappa*dt/h;
t = 0;

% start of integration
resK = zeros(size(gap)); resNa = resK; resCl = resK; resCurrent = 0;
for j=1:ceil(Tmax/ddt)
    current = 0;
    for k=1:round(ddt/dt)
        % K-influx:
        cK(gap) = cK(gap) + kappa*V0;
        
        % Diffusion
        cKTilde = rho2.*cK;
        cNaTilde = rho2.*cNa;    
        cClTilde = rho2.*cCl;    
        cKTilde = cKTilde + DifK*(L*cKTilde-rho2L.*cKTilde);
        cNaTilde = cNaTilde + DifNa*(L*cNaTilde-rho2L.*cNaTilde);
        cClTilde = cClTilde + DifCl*(L*cClTilde-rho2L.*cClTilde);
        cK = cKTilde./rho2;
        cNa = cNaTilde./rho2;
        cCl = cClTilde./rho2;
        current = current + DifK*(sum(rho(outer).*cK(outer))+sum(rho(top).*cK(top)));
        current = current + DifNa*(sum(rho(outer).*cNa(outer))+sum(rho(top).*cNa(top)));
        current = current - DifCl*(sum(rho(outer).*cCl(outer))+sum(rho(top).*cCl(top)));
        
        % Charge relaxation
%         current = current + h/delta*sum(rho(gap).*(cK(gap)+cNa(gap)-cCl(gap)));
%         current = current + sum(rho(vol).*(cK(vol)+cNa(vol)-cCl(vol)));        
        tmp = [cK cNa cCl]*M;
        cK = tmp(:,1);
        cNa = tmp(:,2);
        cCl = tmp(:,3);
%         current = current - h/delta*sum(rho(gap).*(cK(gap)+cNa(gap)-cCl(gap)));
%         current = current - sum(rho(vol).*(cK(vol)+cNa(vol)-cCl(vol)));        
        
    end    
    t(end+1) = t(end) + ddt; 
    
    subplot(311); plot(1e4*rho(bottom), 1e12*cK(bottom));     
    title(['time = ' num2str(1e3*t(end)) ' ms']);     
    ylabel('K^+ [nM]')
    subplot(312); plot(1e4*rho(bottom), 1e12*cCl(bottom)); 
    ylabel('Cl^- [nM]')
    subplot(313); plot(1e4*rho(bottom), 1e12*cNa(bottom)); 
    ylabel('Na^+ [pM]')    
    xlabel('radius [\mum]');
    drawnow; 
    
    resK(:,end+1) = cK(gap);
    resNa(:,end+1) = cNa(gap);
    resCl(:,end+1) = cCl(gap);
    resCurrent(end+1) = 2*pi*delta^2*current;
end
vol = reshape(vol,jump,length(vol)/jump);

% surf(rho(vol),z(vol),cK(vol)*1e12); 
% surf(rho(gap),dt*(1:size(resK,2)),resK'*1e12)
% plot(dt*(1:size(resK,2))*1e3,2*pi*delta*h*rho(gap)'*resK*AvogadroConstant); 
% xlabel('time [ms]'); ylabel('quantity [molecules]')