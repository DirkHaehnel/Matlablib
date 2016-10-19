% program for the integration of the diffusion equation within a cylindrical 
% potential in cylindrical coordinates

% in matrices always:
% rho along first index
% z along second index

nphi = 0:0; % angular components to be considered
nz = 0:20; % z components to be considered
nrho = 1000; 

delta = 0.5e-8; % spatial step size in cm
rho0 = 5*delta;
rhov = rho0 - delta + (0.5:nrho+2)*delta;

DiffConstant = 1e-5; 
BoltzmannConstant = 1.380650324e-23*1e7;
AvogadroConstant = 6.0221419947e23;
ElectronCharge = 1.60217646263e-19*2997924580;
Temperature = 273.15 + 36;
DielectricConstant = 80;
ql = 10*1e7*2*pi*ElectronCharge/DielectricConstant; % statC/cm
ql = ql*2*ElectronCharge/BoltzmannConstant/Temperature;
IonStrength = 1e1*0.5*AvogadroConstant*1e-3; % molecules/ml
kappa = sqrt(8*pi*ElectronCharge^2/DielectricConstant/BoltzmannConstant/Temperature*IonStrength);
psi = real(i*ql*besselh(0,2,-i*kappa*rhov));
dpsi = -real(ql*besselh(1,2,-i*kappa*rhov));
zeta = psi(1);
lam0 = exp(-psi/2).*sqrt(rhov);
lam = (lam0(1:end-2) + lam0(3:end))./lam0(2:end-1) - 2;

M = spdiags([ones(nn,1) -2*ones(nn,1)-lam ones(nn,1)], -1:1, nrho, nrho);
M(
for jz=0:nz
    if jrho==nrho
        M(pos,pos) = lam(jrho) + sqrt(rhov(end))/sqrt(rhov(end-1))*exp((psi(end-1)-psi(end))/2);
        M(pos,nrm) = 1;
    elseif jrho==1
        M(pos,pos) = lam(jrho) + sqrt(rhov(1))/sqrt(rhov(2))*exp((psi(2)-psi(1))/2);;
        M(pos,nrp) = 1;
    end            
end

lam0([1 end]) = [];
psi([1 end]) = [];
rhov([1 end]) = [];

jSource = sum(rhov<5e-8);
% jSink = sum(rhov<10*delta);
% for j=jSink:jSink+5
%     M(j,j) = M(j,j)-1;
% end
conc = zeros(nrho,nz);
conc(jSource,1) = -1;
conc = reshape(conc,prod(size(conc)),1);
conc = M\conc;
conc = reshape(conc,nrho,nz);
conc = conc./((sqrt(rhov).*exp(psi/2))'*ones(1,nz));
% conc = conc/(rhov*conc(:,end)*2*pi*DiffConstant)*1e3; % molecules/cm^3

mim(log10(conc))
% plot(rhov,conc); figure(gcf)

if 0
    k = 20;
    pcolor(zv*1e7,[-flipud(rhov(1:k)'-delta/2); rhov(1:k)'-delta/2]*1e7,log10([flipud(conc0(1:k,:));conc1(1:k,:)]));
    shading interp; 
    axis image
    times16
    patch([0 max(zv) max(zv) 0]*1e7,[-1 -1 1 1]*rho0*1e7, 'b')
    set(gca,'yticklabel',num2str(abs(str2num(get(gca,'yticklabel')))))
    xlabel('\itz\rm [nm]');
    ylabel('\rho [nm]');
    if 1
        h = colorbar('h')
        times16
        set(get(h,'xlabel'),'string','lg(c)','FontName','Times','FontSize',16)    
    else
        h = colorbar
        times16
        set(get(h,'ylabel'),'string','lg(c)','FontName','Times','FontSize',16)    
    end
end

if 0
    pcolor((delta/2+zv(1:end-1))*1e7,[-flipud(rhov(1:20)'); rhov(1:20)']*1e7,log10(-[flipud(diff(conc0(1:20,:)')');diff(conc1(1:20,:)')']/delta));
    shading interp; 
    axis image
    times16
    patch([0 max(zv)-delta/2 max(zv)-delta/2 0]*1e7,[-1 -1 1 1]*rho0*1e7, 'b')
    set(gca,'yticklabel',num2str(abs(str2num(get(gca,'yticklabel')))))
    xlabel('\itz\rm [nm]');
    ylabel('\rho [nm]');
    h = colorbar('h')
    times16
    set(get(h,'xlabel'),'string','lg(-\partialc/\partial\itz\rm) [cm^{-1}]','FontName','Times','FontSize',16)
    
    figure
    semilogy(1e7*rhov,conc0(:,end)/delta,'b',1e7*rhov,conc1(:,end)/delta,'r','linewidth',2)
    times16
    xlabel('\rho [nm]');
    ylabel('-\partialc/\partial\itz\rm [cm^{-1}]')
    axis tight
end


% flux calculation
% plot(-2*pi*delta^2*rhov*diff(conc')'/(rhov(j)*2*pi*delta^2))

