% FCSModel

clear all

radv = (1:0.1:2)*1e3;
cov = 0;
abev = [];
[x,y] = meshgrid(-2:0.05:2,-2:0.05:2);
rr = sqrt(x.^2+y.^2);
pp = angle(x+i*y);
rhofield = [0 3];
zfield = [-4 4];
NAv = 1:0.1:1.3;
n0 = 1.333;
n = 1.333;
n1 = 1.333;
d0 = [];
d = 0;
d1 = [];
lamex = 0.64;
dfoc = 0.2;
focpos = [[0 dfoc 0 0 0];[0 -dfoc 0 0 0]];
pow = 1;
lamem = 0.67;
mag = 60;
av = 100;
zpin = 0e3;
atf = [];
kappa = 1;
lt = [];
fd = 3e3;
resolution = [20 lamex/0.2];
ring = [];
maxm = 10;
pulse = [0.05 12.5]/2;
triplet = 0;
sat = 0;

if 1
    for jrad = 1:length(radv)
        for jNA = 1:length(NAv);
            [jrad jNA]
            NA = NAv(jNA);
            over = [0 radv(jrad) 1];
            [modres(:,:,jrad,jNA), autotime, exc, mdf] = DICFCS(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, ...
                pow, lamem, mag, av, focpos, zpin, atf, kappa, lt, pulse, sat, triplet, resolution, ring, maxm);

            p(:,jrad,jNA) = simplex('GaussFcs',[0.4 0.1],[0 0],[],[],[],av/mag,[lamex lamem]/n,2*dfoc,...
                autotime/1e6,modres(:,1:2,jrad,jNA),modres(:,3,jrad,jNA)/2);

            [err pd(jrad,jNA) c] = GaussFcs(p(:,jrad,jNA),av/mag,[lamex lamem]/n,2*dfoc,...
                autotime/1e6,modres(:,1:2,jrad,jNA),modres(:,3,jrad,jNA)/2);

            vg(jrad,jNA) = GaussDetectionVolume(p(:,jrad,jNA), av/mag, [lamex lamem]/n);
            v1(jrad,jNA) = DetectionVolume(exc.rho,exc.z,mdf.volx1,mdf.voly1);
            v2(jrad,jNA) = DetectionVolume(exc.rho,exc.z,mdf.volx2,mdf.voly2);            
        end
    end
    save FCSModel
end
