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
FaradayConstant = AvogadroConstant*ElectronCharge; % statC/mol
Temperature = 273.15 + 36;
DielectricConstantH2O = 80;
DielectricConstantSiO2 = 3.9;
fac2 = ElectronCharge/BoltzmannConstant/Temperature/2;

% BaSrTiO3 (k = 300)
% Ta2O5 (k = 26); Bandgap = 4.5-eV; Crystalline Structure = Orthorhombic
% ZrO2 (k = 25); Bandgap = 7.8-eV; Crystalline Structure = Monoclinic
% HfO2 (k = 24.5); Bandgap = 5.7-eV; Crystalline Structure = Monoclinic
% Al2O3 (k = 9); Bandgap = 8.7-eV; Amorphous
% Si3N4 (k = 7); Bandgap = 5.1-eV; Amorphous
% SiO2 (Reference) (k = 3.9); Bandgap = 8.9-eV; Amorphous

Tmax = 1; % max integration time
rad = 10e-4; % radius of cell
cyl_thickness = 5e-4; % thickness of enclosing cylinder
cyl_height = 5e-4; % height of enclosing cylinder
cyl_depth = 5e-4; % depth of supporting cylinder (silicon oxide dielectric)
h = 5e-5; % height of gap = FE cell size for integrating the Poisson equation
delta = 5e-5; % FE cell size for integrating the Nernst-Planck equations, should be an odd multiple of h

% diffusion coefficients
DifK = 1.96e-5; 
DifNa = 1.33e-5; 
DifCl = 2.02e-5; 

% bulk concentrations in mol/cm^3
cK0 = 5e-3*1e-3; 
cNa0 = 140e-3*1e-3; 
cCl0 = 145e-3*1e-3; 

% inverse membrane conductivity
kappa = 1/1e10*Coulomb/Volt/FaradayConstant/(4*pi*rad^2); % statC/statV/cm^2/s

% internal potential
V0 = 50e-3*Volt; % statV

jump = round((rad+cyl_thickness)/delta);
[z,rho] = meshgrid(-cyl_depth+delta/2:delta:cyl_height,delta/2:delta:rad+cyl_thickness);
ratio = round(delta/h);
rho = rho(:);
z = z(:);
nn = length(rho);
ind = (1:nn)';
gap = ind(z>0 & z<delta & rho<rad);
bottom = ind(z>0 & z<delta);
outer = ind(rho>(rad+cyl_thickness-delta));
inner = ind(rho>rad & rho<rad+delta & z>delta);
top = ind(z>cyl_height-delta);
down = ind(z<-cyl_depth+delta);
nl = ind-1;
nr = ind+1;
nu = ind+jump;
nd = ind-jump;
nl(rho<delta) = nl(rho<delta)+1;
nr(outer) = ind(outer);
nd(down) = ind(down);
nu(top) = ind(top);
vol = ind(z>0 & rho>rad); 
sil = ind(z<0);
rho2 = sqrt(rho);
sol = [gap;vol];
notsol = setdiff(ind,sol);

% Laplace operator for Poisson integrator
M = spdiags(-4*ones(nn,1),0,nn,nn);
for j=1:nn
    M(j,nr(j)) = M(j,nr(j))+1;
    M(j,nl(j)) = M(j,nl(j))+1;
    M(j,nu(j)) = M(j,nu(j))+1;
    M(j,nd(j)) = M(j,nd(j))+1;
end
L = M;
rho2M = M*rho2;
for j=1:length(gap)
    M(gap(j),nu(gap(j))) = 0;
    M(nu(gap(j)),gap(j)) = 0;
end
for j=1:length(inner)
    M(inner(j),nl(inner(j))) = 0;
    M(nl(inner(j)),inner(j)) = 0;
end
M = M - spdiags(rho2M./rho2,0,nn,nn);

% Laplace operator for Nernst-Planck integrator
for j=1:length(top)
    L(top(j),top(j)) = L(top(j),top(j)) - 1;
end
for j=1:length(outer)
    L(outer(j),outer(j)) = L(outer(j),outer(j)) - 1;
end
for j=1:length(gap)
    L(nu(gap(j)),gap(j)) = 0;
    L(gap(j),nu(gap(j))) = 0;
    L(gap(j),gap(j)) = L(gap(j),gap(j)) + 1;
end
for j=1:length(bottom)
    L(nd(bottom(j)),bottom(j)) = 0;
    L(bottom(j),nd(bottom(j))) = 0;
    L(bottom(j),bottom(j)) = L(bottom(j),bottom(j)) + 1;
end
for j=1:length(inner)
    L(nl(inner(j)),inner(j)) = 0;    
    L(inner(j),nl(inner(j))) = 0;
    L(inner(j),inner(j)) = L(inner(j),inner(j)) + 1;
end
%L(rho<delta & z>0 & z<delta,rho<delta & z>0 & z<delta) = L(rho<delta & z>0 & z<delta,rho<delta & z>0 & z<delta) - 0.5;
L(rho>rad & rho<rad+delta & z>0 & z<delta, rho>rad & rho<rad+delta & z>0 & z<delta) = -2 - h/delta;
L(rho>rad & rho<rad+delta & z>0 & z<delta, rho>rad-delta & rho<rad & z>0 & z<delta) = h/delta;
rho2L = L*rho2;

diel = ones(size(z)); 
diel(z>0) = DielectricConstantH2O; 
diel(z<0) = DielectricConstantSiO2;

% initial conditions
cK = zeros(nn,1); 
cNa = zeros(nn,1); 
cCl = zeros(nn,1); 
qK = zeros(size(gap));
psi = zeros(nn,1);
netCharge = zeros(nn,1);

% relaxation matrix
fac = 4*pi*AvogadroConstant*ElectronCharge^2/BoltzmannConstant/Temperature/DielectricConstantH2O;
cK0 = fac*cK0;
cNa0 = fac*cNa0;
cCl0 = fac*cCl0;

% rescaling and miscellaneous variables
dt = min(0.0002*delta^2./[DifK DifNa DifCl]); % time step in sec
kappa = kappa*dt/h;
t = 0;

% start of integration
dstep = 0.01; % max increment per integration step
respsi = []; resK = []; resNa = []; resCl = [];
for j=1:ceil(Tmax/dt)
    % K-influx:
    cK(gap) = cK(gap) + kappa*(V0-psi(gap)).*(V0>psi(gap));
    
%     ddt = 0;
%     while ddt<dt % integration with self-regulating time step
        psi = 0*psi;
        ee = exp(-fac2.*psi);
        ee2 = ee.^2;
        fm = rho2.*ee;
        Wp = (L*fm)./fm/delta^2;
        fp = rho2./ee;
        Wm = (L*fp)./fp/delta^2;
        cKTilde = fp.*cK;
        cNaTilde = fp.*cNa;    
        cClTilde = fm.*cCl;    
        
        % Nernst-Planck
        tmp = ind(netCharge>0);
        for k=1:length(tmp)
            RM = [-DifK*[cK0+Wp(tmp(k)) cK0 -cK0./ee2(tmp(k))]; ...
                -DifNa*[cNa0 cNa0+Wp(tmp(k)) -cNa0./ee2(tmp(k))]; ...
                DifCl*[cCl0.*ee2(tmp(k)) cCl0.*ee2(tmp(k)) -cCl0-Wm(tmp(k))]];
            conc = expm(dt*RM)*[cKTilde(tmp(k)); cNaTilde(tmp(k)); cClTilde(tmp(k))];
            cKTilde(tmp(k)) = conc(1);
            cNaTilde(tmp(k)) = conc(2); 
            cClTilde(tmp(k)) = conc(3);
            expm(RM)*conc
        end
        cKTilde = cKTilde + dt/delta^2*DifK*L*cKTilde;
        cNaTilde = cNaTilde + dt/delta^2*DifNa*L*cNaTilde;
        cClTilde = cClTilde + dt/delta^2*DifCl*L*cClTilde;        
            
        
        %             cKstep = DifK*(L*cKTilde-Wp.*cKTilde+cK0*fp.*V);
        %             cNastep = DifNa*(L*cNaTilde-Wp.*cNaTilde+cNa0*fp.*V);
        %             cClstep = DifCl*(L*cClTilde-Wm.*cClTilde-cCl0*fm.*V);
        %
        %             tmp = max(-([cKstep(~cKTilde==0)./cKTilde(~cKTilde==0); cNastep(~cNaTilde==0)./cNaTilde(~cNaTilde==0); cClstep(~cClTilde==0)./cClTilde(~cClTilde==0)]));
        %             tmp = dstep/tmp(~isnan(tmp));

        %         tmp = cKstep + cNastep - cClstep;
        %         tmp = abs(dstep*(cKTilde(ind)+cNaTilde(ind)+cClTilde(ind))./tmp(ind));
        %         tmp = min([tmp dt-ddt]);
        %         ddt = ddt + tmp
        %
        %         cKTilde = cKTilde + tmp*cKstep;
        %         cNaTilde = cNaTilde + tmp*cNastep;
        %         cClTilde = cClTilde + tmp*cClstep;

        cK = cKTilde./fp;
        cNa = cNaTilde./fp;
        cCl = cClTilde./fm;
        
        cK(notsol) = 0;
        cNa(notsol) = 0;
        cCl(notsol) = 0;

%         cK(abs(cK)<1e-5*max(abs(cK))) = 0;
%         cNa(abs(cNa)<1e-5*max(abs(cNa))) = 0;
%         cCl(abs(cCl)<1e-5*max(abs(cCl))) = 0;

        % solving the Poisson equation
        netCharge = cK + cNa - cCl;
%        netCharge(gap) = netCharge(gap)*h/delta;
        psi = M\(-delta^2*FaradayConstant*4*pi*rho2.*netCharge./diel);
        psi = psi./rho2;
        clf; 
        subplot(121);
        surf(reshape(rho,jump,length(z)/jump),reshape(z,jump,length(z)/jump),reshape(netCharge,jump,length(z)/jump));
        subplot(122)
        surf(reshape(rho,jump,length(z)/jump),reshape(z,jump,length(z)/jump),reshape(psi,jump,length(z)/jump)); 
        pause
%     end
    t = t + dt; 
    
    subplot(311); plot(1e4*rho(bottom), 1e6*psi(bottom)/Volt); 
    title(['time = ' num2str(1e3*t) ' ms']); 
    ylabel('potential [\muV]')
    subplot(312); plot(1e4*rho(bottom), 1e9*cK(bottom));     
    ylabel('K^+ [\muM]')
    subplot(313); plot(1e4*rho(bottom), 1e9*cNa(bottom), 1e4*rho(bottom), 1e9*cCl(bottom)); 
    ylabel('Na^+/Cl^- [\muM]')    
    xlabel('radius [\mum]');
    drawnow; pause
    
    
    respsi = [respsi psi(gap)];
    resK = [resK cK(gap)];
    resNa = [resNa cNa(gap)];
    resCl = [resCl cCl(gap)];

end

return

mpc(reshape(rho,jump,length(z)/jump),reshape(z,jump,length(z)/jump),reshape(cNa,jump,length(z)/jump))

% surf(reshape(rho,jump,length(rho)/jump),reshape(z,jump,length(rho)/jump),reshape(V,jump,length(rho)/jump)); 
% surf(rho(gap),ddt*(1:size(resK,2)),resK'*1e12)
% plot(ddt*(1:size(resK,2))*1e3,2*pi*delta*h*rho(gap)'*resK*AvogadroConstant); 
% xlabel('time [ms]'); ylabel('quantity [molecules]')

