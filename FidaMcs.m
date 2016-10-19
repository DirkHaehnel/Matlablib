% Programm zur Simulation eines FIDA-Experiments

%clear all
close all

kmax = 1e5;
lamex = 0.64; % excitation wavelength
lamem = 0.67; % emission wavelength
w0 = 0.35; % beam waist 
a0 = 0.15; % confocal aperture parameter
a = 150/60; % confocal aperture radius
d = 0.2; % falf distance between foci
dt = 1e-5; % time interval in microseconds
dif = 10; % diffusion constant in mum^2/sec unit
velo = 500; % flow speed
brightness = 1;
background = 0;

z0 = 5;
N = 2^17; % max number of time channels
L = 1:N;
L2 = L.^2;
stepz = 0.01;
stepx = 0.001;
kzmax = z0/stepz;
kymax = round(5*LaserBeamFit(w0,lamex,0)/stepx);
cnt = 1;
fluo1 = zeros(1,kymax*kzmax);
fluo2 = fluo1;
trans = fluo1;
for kz = 0:kzmax
    scale = LaserBeamFit(w0,lamex,kz*stepz)/LaserBeamFit(w0,lamex,0);
    ds = scale*dt;
    velot = ds*velo;
    xi = sqrt(2*dif*ds);
    z = kz*stepz + cumsum([0 xi*randn(1,N-1)]);
    ymax = kymax*stepx*scale;
    for ky = 0:kymax
        brightness = randi(2);
        NN = ceil(2*sqrt(ymax^2-(ky*stepx*scale)^2)/velot);
        if NN>0
            y = ky*stepx*scale + cumsum([0 xi*randn(1,NN-1)]);
            x = -sqrt(ymax^2-y(1)^2) + velot*(0:NN-1) + cumsum([0 xi*randn(1,NN-1)]);
            %tst(cnt,:) = [x(1) y(1) z(1)];
            
            % fluorescence
            w2 = LaserBeamFit(w0,lamex,z(1:NN)).^2;
            tmp = D1CEF(a0,a,lamem,z(1:NN))./w2*w0^2*scale*brightness;
            f1 = tmp.*exp(-2*(x.^2+(y-d).^2)./w2) + ds*background;
            f2 = tmp.*exp(-2*(x.^2+(y+d).^2)./w2) + ds*background;
            %subplot(121); plot3(x,y,z); axis([-x0 x0 -y0 y0 -z0 z0]); subplot(122); plot(L,f1,L,f2); axis([0 N 0 1]); drawnow
            
            fluo1(cnt) = sum(f1);
            fluo2(cnt) = sum(f2);
            trans(cnt) = ds*sqrt(sum((f1+f2).*L2(1:NN))/sum(f1+f2)-(sum((f1+f2).*L(1:NN))/sum(f1+f2))^2);
            
            if mod(cnt,1e4)==0
                plot([fluo1 fluo2],[fluo2 fluo1],'o');
                axis image;
                title(mnum2str(kz/kzmax,1,4))
                drawnow
            end
            cnt = cnt + 1;
        end
    end
end

[v,bin] = mHist(trans);
[~,pos] = max(v);
bin = bin(pos);
factrans = 0.5;
facint = 0.1;
ind = bin*factrans<trans & trans<bin/factrans & (fluo1-fluo2)./(fluo1+fluo2)*2<facint;
plot([fluo1 fluo2],[fluo2 fluo1],'o',[fluo1(ind) fluo2(ind)],[fluo2(ind) fluo1(ind)],'o'); axis image

figure

zv = 0:0.1:10;
tst = interp1(LaserBeamFit(w0,lamex,zv)/2/velo,zv,trans,'cubic',NaN);
ind = isnan(tst);
fluo1(ind) = [];
fluo2(ind) = [];
trans(ind) = [];
tst(ind) = [];
xv = 0:w0/10:3*w0;
ww = 2*trans*velo;
tmp = D1CEF(a0,a,lamem,tst)./ww*w0;
x = ww.^2/8/d.*log(fluo1./fluo2);
mm = [tmp.*exp(-2*(x-d).^2./ww.^2); tmp.*exp(-2*(x+d).^2./ww.^2)];
ff = sum([fluo1; fluo2].*mm)./sum(mm.^2);
mHist(ff((fluo1+fluo2)>0.1*max(fluo1+fluo2)),1:500)