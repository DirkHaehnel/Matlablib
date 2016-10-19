% AlternativeFCS

close all
% delta = 0.3;
w = 0.3; % lateral extension of 3D-Gaussian volume
pinhole = 50/60; % pinhole radius in object space
lamem = 0.67;
lamex = 0.63;
drho = w/50; dz = drho;
rho = 0:drho:10*w;
z = (0.5:7*w/dz)*dz;
[zz,rr] = meshgrid(z,rho);
cf = D3CEF(rr/lamem*2*pi,zz/lamem*2*pi,pinhole/lamem*2*pi,1.2);
ex = GaussianLaser(rr/lamem*2*pi,zz/lamem*2*pi,w/lamem*2*pi);
volx = cf.*ex;    

ww=2*sqrt(rho.^2*volx./(sum(volx)-volx(1,:)/2));
intens = volx(1,:);

Nsub = 3;
lam = 2^(1/Nsub);
timewin = 1e6;
Ncasc = ceil(log2(timewin));
autotime = lam.^(0:Ncasc*Nsub-1)';
dif = 5e-5; % diffusion constant in mum^2/time unit
modres = autotime;
row = ones(size(z));

% integrating the diffusion equation:
for k = 1:length(autotime)
    tmp = exp(-(row'*z-z'*row).^2/4/dif/autotime(k));
    tmp = tmp + exp(-(row'*z+z'*row).^2/4/dif/autotime(k));        
    tmp = tmp.*(ww'.^2*ww.^2).*(intens'*intens)./(8*dif*autotime(k)+row'*ww.^2+ww'.^2*row)/sqrt(dif*autotime(k));
    modres(k,1) = sum(sum(tmp));
    tmp = tmp.*exp(-2*delta^2./(8*dif*autotime(k)+row'*ww.^2+ww'.^2*row));    
    modres(k,2) = sum(sum(tmp));
    semilogx(autotime(1:k), modres(1:k,:)/modres(1,1)); axis([autotime([1 end])' 0 1.5]); drawnow;
end
modres = modres./modres(1,1);

p = [p simplex('MultiFocusCrossFit',p(:,end),[],[],[],[],autotime, modres(:,1), modres(:,1), modres(:,2))];
p(:,end) = simplex('MultiFocusCrossFit',p(:,end),[],[],[],[],autotime, modres(:,1), modres(:,1), modres(:,2));
% delta^2/4/p(4)