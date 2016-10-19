clear all
close all

% rate_fluo = q/tau0;
% rate_nr = (1-q)/tau0;

tau0 = 3.53; % free lifetime in water
qy = 0.82; % quantum yield
nglass = 1.523; % glass ref index
d0 = 15; % gold thickness
zv = 1:130; % dipole position 
n = 1.3384; % water ref index
d = max(zv);
n1 = n;
d1 = [];
rad = 12.5; % radius of tubulus 
transmission = [516 556]; % transmission band of emission filter FF01-536/40

theta = (0.5:1e3)'/1e3*asin(1.2/nglass);
st = sin(theta);

load Alexa488em.mat
emission = emission(wavelength>transmission(1) & wavelength<transmission(2));
wl = wavelength(wavelength>transmission(1) & wavelength<transmission(2));
gold=load('c:\Joerg\Doc\Nanophotonics&Optics\BerndtMichael\au_2_n-k.txt');
nau = interp1(gold(:,1),gold(:,2)+i*gold(:,3),wl);

exval = [8, 0.91; 27, 1.87];

% lifetime modelling gold
for jw = 1:length(wl)
    lamem = wl(jw);
    n0 = [nglass nau(jw)];

    [lvd,lvu,lpd,lpu,qvd,qvu,qpd,qpu] = LifetimeL(zv/lamem*2*pi,n0,n,n1,d0/lamem*2*pi,d/lamem*2*pi,d1/lamem*2*pi);

    lv(:,jw) = 4/3*n./(qvd+qvu)';
    lp(:,jw) = 4/3*n./(qpd+qpu)';
    

    % The following is unnecessary beacuse for SUB-CRITICAL fluorescence
    % detection, the amount of emitted light into the glass is independent
    % on distance
    %     for jz = 1:length(zv)
    %         [v,pc,ps] = DipoleL(theta,zv(jz)/lamem*2*pi,n0,n,n1,d0/lamem*2*pi,d/lamem*2*pi,d1/lamem*2*pi);
    %         ev(jz,jw) = st'*abs(v).^2;
    %         ep(jz,jw) = st'*(abs(pc).^2+abs(ps).^2)/2;
    %     end
    jw
end

% lifetime modelling glass
[lvd,lvu,lpd,lpu,qvd,qvu,qpd,qpu] = LifetimeL(zv/lamem*2*pi,nglass,n,n1,[],d/lamem*2*pi,d1/lamem*2*pi);

gp = (qy./lp*emission/sum(emission) + 1-qy)/tau0;
gv = (qy./lv*emission/sum(emission) + 1-qy)/tau0;
phi = (0.5:360)/180*pi;
ind = zv<max(zv)-2*rad;
zz = zv(ind);
for j=1:length(zz)
    zi = zv(j) + rad + rad*cos(phi);
    gi = interp1(zv,(gv+2*gp)/3,zi);
    gg(j) = sum(1./gi.^2)/sum(1./gi);
end

wv = tau0./(qy./(4/3*n./(qvd+qvu)') + 1-qy);
wp = tau0./(qy./(4/3*n./(qpd+qpu)') + 1-qy);
ww = tau0./(qy./(4*n./(2*qpd+2*qpu+qvd+qvu)') + 1-qy);

plot(zz,ww(ind),zz,gg,zz,ones(size(zz))*tau0,'--',exval(:,1),exval(:,2),'sk')
grid
axis([0 max(zz) 0 4])
xlabel('distance (nm)')
ylabel('lifetime (ns)')

save 'FlicFlim2009-10-04'