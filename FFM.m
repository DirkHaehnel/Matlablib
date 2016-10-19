nmax = 100; % number of emitters
rad = 50; %radius of circle
br = 5; %number of photons per emitter per time step
kp = 0.01; % forward rate constant into dark state
km = 0.001; % backward rate constant into bright state
tmax = 1e4; % number of time steps

NA = 1.2;
al = pi/2; be = 0;
mag = 1000;
lamem = 0.58;
nh2o = 1.33;
focpos = 0;
[intx, inty, intz, rho, phi, fxx0, fxx2, fxz, byx0, byx2, byz] = ...
    SEPDipole([0 2], 0, NA, nh2o, nh2o, nh2o, [], 0, [], lamem, mag, focpos);

int = SEPImage(pi/2,0,305,1,rho,fxx0,fxx2,fxz,byx0,byx2,byz);
int = int + SEPImage(pi/2,pi/2,305,1,rho,fxx0,fxx2,fxz,byx0,byx2,byz);
int = int + SEPImage(0,0,305,1,rho,fxx0,fxx2,fxz,byx0,byx2,byz);
[mx, my, wx, wy, om, amp] = GaussEllipse([],[],int);

clear intx inty intz rho phi fxx0 fxx2 fxz byx0 byx2 byz

[xx,yy] = meshgrid(-304.5:305,-304.5:305);
phi = 2*pi*rand(nmax);
px = rad*cos(phi); py = rad*sin(phi);

for j=1:nmax
    tmp = exp(-((xx-px(j)).^2 + (yy-py(j)).^2)/2/wx.^2);
    for kx=1:61
        for ky=1:61
            im0(kx,ky,j) = sum(sum(tmp((kx-1)*10+1:kx*10,(ky-1)*10+1:ky*10)));
        end
    end
    im0(:,:,j) = br*im0(:,:,j)/sum(sum(im0(:,:,j)));
end

state_vec = ones(nmax,1); % state vector of emitters
for k=1:tmax
    state_vec = state_vec - (state_vec==1).*(rand(nmax,1)<=kp) + (state_vec==0).*(rand(nmax,1)<=km);
    im(:,:,k) = poissrnd(im0(:,:,1));
    for j=2:nmax
        im(:,:,k) = im(:,:,k) + poissrnd(im0(:,:,j));
    end
end

