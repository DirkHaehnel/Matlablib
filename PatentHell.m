% PatentHell

NA = 1.4;
n = 1.51;
over = [NA*3e3 inf 5e3];
lamex = 0.514;
[fx0, fx2, fz, rho, z, exc, excx, excz] = NewExc([-0.00125 1.5], 0, NA, n, n, n, [], 0, [], lamex, over, 0, [], lamex/0.0025);
exc = exc(:,1)';
rho = rho(:,1)';
clear fx0 fx2 fz z excx ecxz


% Acridin Orange
lambda = 470e-9;
extinction = 5e4;
tau = 2.7e-9;
kisc = 0.1/tau;
kph = 1/3e-6;

NatConst;
hnu = PlanckConstant*SpeedOfLight/lambda;

power = 5e-5; % Anregungsleistung in W/cm^2
area = pi*(1e-4*0.18)^2; 
intensity = 2*power/area;

sigma = extinction*log(10)*1e3/AvogadroConstant;

exc0 = sigma*intensity/hnu;
tv = [0.05 0.1 0.2 0.4 0.8 1.6]*1e-6;

% Fig.1a
z = exc'*exc0;
for j=1:length(tv)
    z(:,j+1) = Hell(exc0*exc,tv(j),tau,kisc,kph);
end
plot(rho,z/exc0,'k');
rhov = [0 0.04 0.08 0.16 0.32];
hold on
for j=2:length(rhov)
    plot(rhov(j)*[1 1],[0 1],':k');
end
hold off

% Fig.1b
t = (0:0.01:3)*1e-6;
x = exc(1);
for j=2:length(rhov)
    x(j) = exc(chop(rho,3)==rhov(j));
end
zz = Hell(exc0*x,t,tau,kisc,kph)';
plot(t,zz./(ones(size(zz,1),1)*max(zz)),'k')


% Fig.2
y = -1.:0.0025:1.;
[x,y] = meshgrid(y,y);
x = sqrt(x.^2+y.^2);
z1 = interp1(rho,exc,x,'cubic');
clear x y

y = -1.:0.0025:1.;
im = zeros(length(y));
im(:,(size(im,2)-1)/2) = 1;
im0 = mconv2(im,z1);
clear imm
tv = [0 0.05 0.1 0.2 0.4 0.8 1.6 3.2]*1e-6;
tauv = [0.75 1.5]*1e-6;

for j=2:length(tv)
    t = [tv(j-1) tv(j)];
    imm(:,:,j-1) = mconv2(im,Hell(exc0*z1,t,tau,kisc,kph));
    MM(j-1,:) = [diff(t) -diff(exp(-t'*(1./tauv))).*tauv];
end
for j=1:size(im,2) tst1(:,j)=lsqnonneg(MM,squeeze(imm(round(end/2),j,:))); end

plot(y,im0((end-1)/2,:)/max(im0(:)),y,tst1(2,:)/max(tst1(2,:)),[-1 1],[1 1]/exp(2),':')
axis([-0.75 0.75 0 1])
xlabel('\itx\rm [\mum]');
ylabel('rel. Amplitude');

sqrt([sum(y.^2.*im0((end-1)/2,:))/sum(im0((end-1)/2,:)),sum(y.^2.*tst1(2,:))/sum(tst1(2,:))])
