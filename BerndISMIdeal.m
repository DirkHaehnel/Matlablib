close all

NA = 1.2;
fd = 3e3;
n0 = 1.33;
n = 1.33;
n1 = 1.33;
d0 = [];
d = [];
d1 = [];
lamex = 0.56;
over = 8e3;
focpos = 0;
defoc = [];
lamem = 0.59;
mag = 60;
zpin = 0;
atf = [];
resolution = 500;
ring =  'besselj(1,rad/5)./(1e-10+rad/5)';
bild = 0;

rhofield = [0 0.5] - lamex/2/resolution;
zfield = [0 2*lamex/resolution];

avv = (0.02:0.04:1.4)*60/2;

[x,y] = meshgrid(0:lamex/resolution:1,0:lamex/resolution:1);
rr = sqrt(x.^2+y.^2);
bead_profile = real(sqrt(0.1^2-rr.^2));

for j=1:length(avv)
    av = avv(j);
    [exc, wave, mdf] = MDFConfocalMicroscopy(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, defoc, av, lamem, mag, zpin, atf, resolution, ring, [], [], [], [], [], bild);
    tmp = interp1(exc.rho(:,1), mdf.volx(:,1,1)+mdf.voly(:,1,1), rr, 'cubic', 0);
    subplot(121)
    [bla, bla, wx1(j), wy1(j)] = Gauss2D(x(1,:),y(:,1),tmp);
    tmp = mConv2(tmp,bead_profile);
    subplot(122)
    [bla, bla, wx2(j), wy2(j)] = Gauss2D(x(1,:),y(:,1),tmp);
    drawnow
end