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
DielectricConstantSiO2 = 80; %3.9;
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
h = 1e-5; % height of gap = FE cell size for integrating the Poisson equation
delta = h; % FE cell size for integrating the Nernst-Planck equations, should be an odd multiple of h

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
V0 = 50e-3*1e6/SoL; % statV

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
nl(rho<h) = nl(rho<delta)+1;
nr(outer) = ind(outer);
nd(down) = ind(down);
nu(top) = ind(top);
vol = ind(z>0 & rho>rad); 
sil = ind(z<0);
rho2 = sqrt(rho);

sol = [gap;vol];
rho_sol = rho(sol);
z_sol = z(sol);
ind_sol = ind(sol);
nn_sol = length(rho_sol);
tmp = z_sol + 2*max(z_sol(:))*rho_sol;
[tmp, tmp, sol2psi] = unique(tmp);
rho_sol = rho_sol(tmp);
z_sol = z_sol(tmp);
rho_sol2 = sqrt(rho_sol);
nn_sol = length(rho_sol);
ind_sol = (1:nn_sol)';
gap_sol = ind_sol(z_sol>0 & z_sol<delta & rho_sol<rad);
bottom_sol = ind_sol(z_sol>0 & z_sol<delta);
outer_sol = ind_sol(rho_sol>rad+cyl_thickness-delta);
inner_sol = ind_sol(z_sol>delta & rho_sol>rad & rho_sol<rad+delta);
top_sol = ind_sol(z_sol>cyl_height-delta);
nl_sol = ind_sol-1;
nr_sol = ind_sol+1;
jump_sol = round(cyl_thickness/delta);
nu_sol(1:round(rad/delta)) = ind_sol(1:round(rad/delta));
nd_sol(1:round((cyl_thickness+rad)/delta)) = ind_sol(1:round((cyl_thickness+rad)/delta));
nu_sol(round(rad/delta)+1:nn_sol) = ind_sol(round(rad/delta)+1:nn_sol) + jump_sol;
nd_sol(round((cyl_thickness+rad)/delta)+1:nn_sol) = ind_sol(round((cyl_thickness+rad)/delta)+1:nn_sol) - jump_sol;
nl_sol(1) = 1;
nr_sol(outer_sol) = ind_sol(outer_sol);
nl_sol(inner_sol) = ind_sol(inner_sol);
nu_sol(top_sol) = ind_sol(top_sol);

% solution-FE to potential-FE mapping matrix
psi2sol = spdiags(zeros(nn_sol,1),0,nn_sol,nn);
tmp = ind(sol);
for j=1:nn_sol
    psi2sol(j,tmp(sol2psi==j)) = 1/length(tmp(sol2psi==j));
end

% Laplace operator for Nernst-Planck integrator
L = spdiags(-4*ones(nn_sol,1),0,nn_sol,nn_sol);
for j=1:nn_sol
    L(j,nr_sol(j)) = L(j,nr_sol(j))+1;
    L(j,nl_sol(j)) = L(j,nl_sol(j))+1;
    L(j,nu_sol(j)) = L(j,nu_sol(j))+1;
    L(j,nd_sol(j)) = L(j,nd_sol(j))+1;
end
L(ind_sol(rho_sol>rad & rho_sol<rad+delta & z_sol<delta),ind_sol(rho_sol>rad & rho_sol<rad+delta & z_sol<delta)) = -2 - h/delta;
L(ind_sol(rho_sol>rad & rho_sol<rad+delta & z_sol<delta),ind_sol(rho_sol>rad-delta & rho_sol<rad  & z_sol<delta)) = h/delta;
rho_sol2L = L*rho_sol2;

% Laplace operator for Poisson integrator
M = spdiags(zeros(nn,1),0,nn,nn);
for j=1:nn
    if z(j)>0
        M(j,j) = -4*DielectricConstantH2O;
        M(j,nr(j)) = M(j,nr(j))+DielectricConstantH2O;
        M(j,nl(j)) = M(j,nl(j))+DielectricConstantH2O;
        M(j,nu(j)) = M(j,nu(j))+DielectricConstantH2O;
        M(j,nd(j)) = M(j,nd(j))+DielectricConstantH2O;
    else
        M(j,j) = -4*DielectricConstantSiO2;
        M(j,nr(j)) = M(j,nr(j))+DielectricConstantSiO2;
        M(j,nl(j)) = M(j,nl(j))+DielectricConstantSiO2;
        M(j,nu(j)) = M(j,nu(j))+DielectricConstantSiO2;
        M(j,nd(j)) = M(j,nd(j))+DielectricConstantSiO2;
    end    
end
rho2M = M*rho2;
% for j=1:length(outer)
%     if z(outer(j))>0
%         M(outer(j),outer(j)) = M(outer(j),outer(j)) - DielectricConstantH2O;
%     else
%         M(outer(j),outer(j)) = M(outer(j),outer(j)) - DielectricConstantSiO2;
%     end
% end
for j=1:length(top)
    M(top(j),top(j)) = M(top(j),top(j)) - DielectricConstantH2O;
end
for j=1:length(down)
    M(down(j),down(j)) = M(down(j),down(j)) - DielectricConstantSiO2;
end
for j=1:length(gap)
    M(gap(j),nu(gap(j))) = 0;
end
for j=1:length(inner)
    M(inner(j),nl(inner(j))) = 0;
end
M = M - spdiags(rho2M./rho2,0,nn,nn);

% initial conditions
cK = zeros(nn_sol,1); 
cNa = zeros(nn_sol,1); 
cCl = zeros(nn_sol,1); 
qK = zeros(size(gap_sol));
psi = zeros(nn_sol,1);
netCharge = zeros(nn,1);

% rescaling and miscellaneous variables
dt = 1e-3*min(0.2*delta^2./[DifK DifNa DifCl]); % time step in sec
DifK = DifK/delta^2; 
DifNa = DifNa/delta^2; 
DifCl = DifCl/delta^2; 
kappa = kappa*dt/h;
t = 0;

% start of integration
dstep = 0.2; % max increment per integration step
respsi = []; resK = []; resNa = []; resCl = [];
for j=1:ceil(Tmax/ddt)
    % K-influx:
    if j==1 & k==1 
        cKTilde(gap_sol) = cKTilde(gap_sol) + kappa*fp(gap_sol).*(V0-psi(gap_sol)).*(V0>psi(gap_sol));
    end
    
    ddt = 0;
    while ddt<dt % integration with self-regulating time step
        fm = rho_sol2.*exp(-fac2.*psi);
        Wp = (L*fm)./fm;
        fp = rho_sol2.*exp(fac2*psi);
        Wm = (L*fp)./fp;
        V = ElectronCharge/Temperature/BoltzmannConstant*(L*(rho_sol2.*psi) - psi.*rho_sol2L)./rho_sol2;
        cKTilde = fp.*cK;
        cNaTilde = fp.*cNa;    
        cClTilde = fm.*cCl;    
        
        % Nernst-Planck
        cKstep = DifK*(L*cKTilde-Wp.*cKTilde+cK0*fp.*V);
        cNastep = DifNa*(L*cNaTilde-Wp.*cNaTilde+cNa0*fp.*V);
        cCLstep = DifCl*(L*cClTilde-Wm.*cClTilde-cCl0*fm.*V);
      
        tmp = min(dstep*[min(cKTilde./cKstep) min(cNaTilde./cNastep) min(cClTilde./cClstep)]);
        tmp = min([tmp dt-ddt]);
        ddt = ddt + tmp
        
        cKTilde = cKTilde + tmp*DifK*(L*cKTilde-Wp.*cKTilde+cK0*fp.*V);
        cNaTilde = cNaTilde + tmp*DifNa*(L*cNaTilde-Wp.*cNaTilde+cNa0*fp.*V);
        cClTilde = cClTilde + tmp*DifCl*(L*cClTilde-Wm.*cClTilde-cCl0*fm.*V);

        cK = cKTilde./fp;
        cNa = cNaTilde./fp;
        cCl = cClTilde./fm;
        
        % boundary conditions
        cK(outer_sol) = 0;
        cNa(outer_sol) = 0;
        cCl(outer_sol) = 0;
        cK(top_sol) = 0;
        cNa(top_sol) = 0;
        cCl(top_sol) = 0;

        % solving the Poisson equation
        netCharge(sol) = cK(sol2psi) + cNa(sol2psi) - cCl(sol2psi);
        psi = M\(-h^2*FaradayConstant*4*pi*rho2.*netCharge);
        psi = psi./rho2;
        psi = psi2sol*psi;
    end
    t = t + dt; 
    
    subplot(311); plot(1e4*rho_sol(bottom_sol), 1e3*psi(bottom_sol)/Volt, 1e4*rho_sol(bottom_sol), 1e3*ones(length(bottom_sol),1)*V0/Volt); 
    title(['time = ' num2str(1e3*t) ' ms']); 
    ylabel('potential [mV]')
    subplot(312); plot(1e4*rho_sol(bottom_sol), 1e12*cK(bottom_sol));     
    ylabel('K^+ [nM]')
    subplot(313); plot(1e4*rho_sol(bottom_sol), 1e15*cNa(bottom_sol), 1e4*rho_sol(bottom_sol), 1e15*cCl(bottom_sol)); 
    ylabel('Na^+/Cl^- [pM]')    
    xlabel('radius [\mum]');
    drawnow; % pause
    
    respsi = [respsi psi(gap_sol)];
    resK = [resK cK(gap_sol)];
    resNa = [resNa cNa(gap_sol)];
    resCl = [resCl cCl(gap_sol)];

end

% surf(rho(vol),z(vol),cK(vol)*1e12); 
% surf(rho(gap),ddt*(1:size(resK,2)),resK'*1e12)
% plot(ddt*(1:size(resK,2))*1e3,2*pi*delta*h*rho(gap)'*resK*AvogadroConstant); 
% xlabel('time [ms]'); ylabel('quantity [molecules]')