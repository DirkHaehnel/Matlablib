% FCSFit

close all
AvogadroConstant = 6.0221419947e23;

path = 'D:\Joerg\Doc\Fcs\2Focus\GDCL\2006-04-11\'
names = dir([path 'Atto655*.mat']);

for j=1:length(names)
    pos = findstr(names(j).name,'P');
    zd(j) = str2num(names(j).name(pos+1:pos+1));
end
[zd ind] = sort(zd);
names = names(ind);

av = 100/58;
lam = [0.64 0.67]/1.33;
dist = 0.4;

global pd
expflag = [];
clear p v dc trip conc TT
for j=1:length(names)
    [auto,t] = FCSCrossRead([path names(j).name],[1e-6 5e-1]);
    eval([path names(j).name ' parameters']); 
    TT(j) = str2num(parameters.Temperature);

    pd = [1 1e-3 1e-3 1e-3];
    p(:,j) = simplex('GaussFCS',[0.37 0.12],zeros(1,2),[],[],[],av,lam,dist,1e2*t,auto(:,1:2),auto(:,3),expflag);

    [err, pd, c] = GaussFcs(p(:,j),av,lam,dist,1e2*t,auto(:,1:2),auto(:,3),expflag);

    v(j) = GaussDetectionVolume(p(:,j),av,lam);
    dc(j) = 1./pd(1)*1e-6
    if ~isempty(expflag)
        trip(:,j) = pd(2:end)*1e-2;
    end
    conc(j,:) = c(2,:)./c(1,:)/1e-15/AvogadroConstant;
    
end

% save GaussFCSFit zd dc p names av lam dist

return

load viscosity
vis = interp1(T,visc(1,:),273.15+TT);
vis0 = interp1(T,visc(1,:),273.15+25);
dcc = dc.*vis./(TT+273.15)*(25+273.15)/vis0


ref = [1.331 1.3692 1.388 1.4057 1.4138 1.4227 1.4299 1.4361 1.4408];
visc_theo = [1 1.11 1.199 1.329 1.409 1.518 1.617 1.713 1.797];
mol_theo =  [0 2.209 3.3189 4.418 4.940 5.513 5.973 6.348 6.632];
mol_ref =  (ref-min(ref))/(max(ref)-min(ref))*max(mol_theo);
visc_ref = interp1(mol_theo,visc_theo,mol_ref,'cubic');


