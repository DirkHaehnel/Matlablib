% temporal deacy at certain points along the axis
close
a = 100e-7; % radius of tube in cm
DiffConstant = 1e-6; 
AvogadroConstant = 6.0221419947e23;
tau = 1e-3; % time step in s
tinflux = 300*tau;
kappa = 5; % GC activity in 1/s
tv = [3*10.^(0:2) 1e3 1e4]*tau;
z = (-2.5e3:1e1:2.5e3)'*1e-6;
col = ones(size(z));
row = ones(size(tv));
conc = col*row*0;
for j=-3:3
    zz = z+2*j*max(z);
    t = tv;
    conc = conc + (col*sqrt(t/pi/DiffConstant)).*exp(-zz.^2*(1./t/DiffConstant/4)) + abs(zz*row/2/DiffConstant).*erf(abs(zz)*sqrt(1./t/DiffConstant/4));
    t = max([eps*ones(size(tv));tv-tinflux]);
    conc = conc - (col*sqrt(t/pi/DiffConstant)).*exp(-zz.^2*(1./t/DiffConstant/4)) - abs(zz*row/2/DiffConstant).*erf(abs(zz)*sqrt(1./t/DiffConstant/4));
end
plot(z*1e4,kappa*conc/AvogadroConstant*1e3/pi/a^2*1e9)
ax = axis;
ax(3) = 0;
axis(ax)
% clear s
% for j=1:length(tv)
%     s{j} = ['\itz\rm = ' num2str(tv(j)) ' s'];
% end
% legend(s)
box off
xlabel('z (\mum)');
ylabel('concentration (nM)')
ax = axis;
% text(ax(1) + 0.1*diff(ax(1:2)),ax(3) + 0.7*diff(ax(3:4)),{'\itD\rm = 10^{-6} cm^2/s', '', '\itk\rm_{GC} = 5/s', '', '\itT\rm_{GC} = 300 ms'})

return

% program for the integration of the diffusion equation in cylindrical
% coordinates

% in matrices always:
% rho along first index
% z along second index

close

a = 100e-7; % radius of tube in cm
DiffConstant = 1e-6; 

BoltzmannConstant = 1.380650324e-23*1e7;
AvogadroConstant = 6.0221419947e23;
ElectronCharge = 1.60217646263e-19*2997924580;
Temperature = 273.15 + 36;
DielectricConstant = 80;

load BesselJPrimeZeros
% nbessel = 120+1; % number of Bessel functions to be considered
% nroots = 100; % number of Bessel function roots to be considered
% for jb=1:nbessel
%     [tmp, tmp] = BesselRoots(jb-1,nroots);
%     for j=1:nroots 
%         tmp(j)=fzero(inline(['besselj(' int2str(jb-2) ',x)-bessel(' int2str(jb) ',x)']),[tmp(j)-1 tmp(j)+1],optimset('TolX',1e-15)); 
%     end    
%     rb(jb,:) = tmp;
% end

% for jb=1:nbessel
%     cnalpha(jb,:) = 1./(pi^2*a^2*(1-(jb-1)^2./(rb(jb,:).^2+(jb==1))).*besselj(jb-1,rb(jb,:)));
% end
% cnalpha(1,:) = cnalpha(1,:)/2;
% cnalpha = pi/2/DiffConst*cnalpha./rb;

cnalpha = rb./4/pi/DiffConstant/a./(rb.^2-(0:nbessel-1)'.^2*ones(1,nroots));
cnalpha(1,:) = cnalpha(1,:)/2;

dphi = 2*pi/(nbessel-1);
phi = (-pi:dphi:pi)';
z = (0:5:5000)*1e-7;
row = ones(size(z));
col = ones(nroots,1);

close
tau = 1e-6; % time step in s
tv = tau*exp((0:100)/100*log(1e5));
xx = ones(size(phi))*z; yy = -a*cos(phi)*ones(size(z)); zz = a*sin(phi)*ones(size(z));
for jt=1:length(tv)
    t = tv(jt);
    tmpz = sqrt(z.^2/4/DiffConstant/t);
    conc = 1/2/pi/a^2/DiffConstant*ones(size(phi))*(2*sqrt(DiffConstant*t/pi)*exp(-tmpz.^2)+z.*(erf(tmpz)-1));
    for jb=1:nbessel
        conc = conc + cos((jb-1)*phi)*sum((cnalpha(jb,:)'*row).*(exp(-rb(jb,:)'*z/a).*(erf(rb(jb,:)'*row*sqrt(DiffConstant*t)-col*tmpz)+1)+...
            exp(rb(jb,:)'*z/a).*erfc(rb(jb,:)'*row*sqrt(DiffConstant*t)+col*tmpz)));
    end
    surf([-fliplr(xx) xx],[yy yy],[zz zz],[fliplr(conc) conc]); axis image; view([-25 20]);
    camlight left 
    camlight headlight
    axis off
    line([0 1e-4],-12*[a a],-0*[a a],'linewidth',3,'color','k')
    text(0,-20*a,-0*a,'1 \mum')
    title({['time = ' mint2str(t*1e6,3,' ') ' \mus'],['max. conc = ' mnum2str(max(conc(:))*1e3/AvogadroConstant*1e9,2,2) ' nM']})
    shading interp
    drawnow
    eval(['print -dpng tmp' mint2str(jt,3)])
end

break

% cross section at z = 0
rho = (0:delta:a)';
col = ones(size(rho));
conc = 1/2/pi^2/a^2*ones(length(phi),length(rho));    
for jb=1:nbessel
    conc = conc + cos((jb-1)*phi)*sum(((col*(cnalpha(jb,:).*exp(-DiffConstant/a^2*rb(jb,:).^2*t))).*besselj(jb-1,rho*rb(jb,:)/a))');
end
surf(cos(phi)*rho',sin(phi)*rho',conc)

return

% temporal deacy at certain points along the circumference
close
tau = 1e-7; % time step in s
tv = (5:1e2)*tau;
for jb=1:nbessel
    M(jb,:) = cnalpha(jb,:).*besselj(jb-1,rb(jb,:));
end
phi = (0:25:50)'*1e-7/a;
clear conc
for jt=1:length(tv)
    t = tv(jt);
    tmp = 1/2/pi^2/a^2*ones(size(phi));    
    tmp = tmp + cos(phi*(0:nbessel-1))*sum((M.*exp(-DiffConstant/a^2*rb.^2*t))')';
    conc(:,jt) = tmp/sqrt(DiffConstant*t/pi);
end
semilogy(tv*1e6,conc/AvogadroConstant*1e3)
axis([0 10 1e-9 1e-2])
for j=1:length(z)
    s{j} = ['\ita\rm\cdot\phi = ' num2str(z(j)*1e7) ' nm'];
end
legend(s)
xlabel('time [\mus]');
ylabel('concentration [M]')


break
% temporal deacy at certain points along the circumference 2
close
tau = 1e-7; % time step in s
tv = 5*exp((0:0.1:1)*log(1e3))*tau;
for jb=1:nbessel
    M(jb,:) = cnalpha(jb,:).*besselj(jb-1,rb(jb,:));
end
dphi = 2*pi/(nbessel-1);
phi = (-pi:dphi:pi)';
clear conc
for jt=1:length(tv)
    t = tv(jt);
    tmp = 1/2/pi^2/a^2*ones(size(phi));    
    tmp = tmp + cos(phi*(0:nbessel-1))*sum((M.*exp(-DiffConstant/a^2*rb.^2*t))')';
    conc(:,jt) = tmp/sqrt(DiffConstant*t/pi);
end
len = size(conc,2);
semilogy(phi*a*1e7,conc(:,1)/AvogadroConstant*1e3,'r')
hold on
for j=2:len
    semilogy(phi*a*1e7,conc(:,j)/AvogadroConstant*1e3,'color',[1-j/len 0 j/len]);
end
axis([-pi*a*1e7 pi*a*1e7 1e-10 1e-2])
for j=1:len
    s{j} = ['\itt\rm = ' mnum2str(tv(j)*1e6,3,1) ' \mus'];
end
legend(s)
xlabel('\ita\rm\cdot\phi [nm]');
ylabel('concentration [M]')

break 
% Rich & Karpen
tau = 1e-7; % time step in s
tv = exp((0:0.001:1)*log(1e7))*tau;
loglog(tv,2./(4*pi*DiffConstant*tv).^(3/2).*exp(-25e-7.^2/4/DiffConstant./tv)/AvogadroConstant*1e3,tv,1/2/pi/DiffConstant/25e-7*(1-erf(25e-7./sqrt(4*DiffConstant*tv)))/AvogadroConstant*1e3)
xlabel('time [s]');
ylabel('concentration [M]')


break
% Surface integral
len = 100e-4; 
load BesselJPrimeZeros
tau = 1e-7; % time step in s
tv = exp((0:0.001:1)*log(1e7))*tau;
semilogx(tv,1/pi/a^2/len+sum(exp(-DiffConstant*rb(1,:)'.^2/a^2*tv)))

mm = (-100:100)';
tst = sum(erf((2*mm+1)*len*sqrt(1./tv)/sqrt(4*DiffConstant)) - erf((2*mm-1)*len*sqrt(1./tv)/sqrt(4*DiffConstant)));

