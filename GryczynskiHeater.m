close all

NatConst
conc = 1e-3; %in M
conc = conc*AvogadroConstant/1e15; 
sigma = 150e3; % extinction in l/M/cm
sigma = sigma*log(10)*1e3/AvogadroConstant*1e8;
qyield = 1; % quantum yield of heat conversion
pow = 0.1; % total power in W
heat_cap = 4.184/1e12; % volumetric heat capacity J/K/mum^3
dt = 1e-9; % time step in s

NA = 1.2;
fd = 3e3;
n0 = 1.333;
n = 1.333;
n1 = 1.333;
d0 = [];
d = [];
d1 = [];
lamex = 0.75;
over = [0 5e3];
focpos = 10;
atf = [];

if ~exist('exc')==1
    resolution = 30;
    rhofield = [0 10];
    zfield = [0 20];
    exc = GaussExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, atf, resolution, [], 1);
    s0 = 2*pi*diff(exc.rho(1:2,1))*max(exc.rho(:,1)'*(abs(exc.fxc(:,:,1)).^2 + abs(exc.fyc(:,:,1)).^2 + abs(exc.fzc(:,:,1)).^2));
end
f1 = (abs(exc.fxc(:,:,1)).^2 + abs(exc.fyc(:,:,1)).^2 + abs(exc.fzc(:,:,1)).^2)/s0; % only rotationally symmetric components
f1 = dt*pow*qyield*sigma*conc*diff(exc.z(1,1:2))*diff(exc.rho(1:2,1))*f1/heat_cap;

% integrating the thermal diffusion equation:

% water heat conductivity = 0.58 W/m/K
dif = 0.58/heat_cap/1e6; % thermal diffusion constant in mum^2/time unit

rhom = exc.rho;
z = exc.z;
dz = z(1,2)-z(1,1);
z = ones(size(z,2),1)*z(1,:);
drho = rhom(2,1)-rhom(1,1); 
rho = rhom(:,1)*ones(1,size(rhom,1));
block = ones(size(rho(:,1)))*rho(:,1)';
ind1 = exc.rho(:,1)<=2; ind2 = abs(exc.z(1,:)-focpos)<=2; % region of interest
plot_time = 0.1e-6*2.^(0:6);
ff1 = 0*f1;
ff2 = 0*f1;
T = 1e4;
for k = 1:2*T
    s4dt = sqrt(4*dif*k*dt);
    tmpz = (erf((z-z'+dz/2)/s4dt) - erf((z-z'-dz/2)/s4dt))/2;
    rr = rho(:,1)*rho(:,1)'/2/dif/k/dt;
    rpr = (rho-rho').^2/4/dif/k/dt;
    tmp = besseli(0, rr, 1).*exp(-rpr);
    tmp(isnan(tmp)) = 0;
    tmp = block.*tmp/2/dif/k/dt*drho;
    ss = sum(tmp')';
    tmp(ss>1,:) = diag(1./ss(ss>1))*tmp(ss>1,:);
	ff1 = ff1 + tmp*f1*tmpz;
    if k>T
        s4dt = sqrt(4*dif*(k-T)*dt);
        tmpz = (erf((z-z'+dz/2)/s4dt) - erf((z-z'-dz/2)/s4dt))/2;
        rr = rho(:,1)*rho(:,1)'/2/dif/(k-T)/dt;
        rpr = (rho-rho').^2/4/dif/(k-T)/dt;
        tmp = besseli(0, rr, 1).*exp(-rpr);
        tmp(isnan(tmp)) = 0;
        tmp = block.*tmp/2/dif/(k-T)/dt*drho;
        ss = sum(tmp')';
        tmp(ss>1,:) = diag(1./ss(ss>1))*tmp(ss>1,:);
        ff2 = ff2 + tmp*f1*tmpz;
    end
    if mod(k,1e1)==0
        surf([-flipud(exc.rho(ind1,ind2)); exc.rho(ind1,ind2)], ...
            [exc.z(ind1,ind2); exc.z(ind1,ind2)], ....
            [flipud(ff1(ind1,ind2)-ff2(ind1,ind2)); ff1(ind1,ind2)-ff2(ind1,ind2)])
        shading flat
        ax = axis;
        axis([ax(1:5) 1.5])
        set(gca,'xticklabel',abs(get(gca,'xtick')))
        xlabel('\rho (\mum)')
        ylabel('\itz\rm (\mum)')
        zlabel('temperature difference (K)');
        title(['time = ' mnum2str(k*dt*1e6) ' \mus'])
        drawnow
        if ~isempty(intersect(k,round(plot_time/dt))) || ~isempty(intersect(k-T,round(plot_time/dt)))
            eval(['print -dpng -r300 -zbuffer GrycziHeater' mint2str(1e9*k*dt,6) 'ns_A']);
        end
    end
    res(k) = max(ff1(:)-ff2(:));
end
% save GrycziHeater res

plot(dt*1e6*(1:length(res)),res)
xlabel('time (\mus)')
ylabel('peak tempeature (K)')

