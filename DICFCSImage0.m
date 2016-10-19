radv = [1.25 1.5 2 3 4]*1e3;
covv = 0:10;
astv = 0:0.1:0.3;

dist = 0.40;
rhofield = [0 1.5];
zfield = [9 23];
NA = 1.14;
fd = 3e3;
n0 = 1.333;
n = n0;
n1 = n0;
d0 = [];
d = 0;
d1 = [];
lamex = 0.64;
focposexc = [15 dist/2 0 0 0];
pow = 1;
lamem = 0.67;
mag = 60;
av = 100;
focposdet = 15;
zpin = 0e3;
kappa = 0;
lt = [];
pulse = [0.05 25]/2; % laser characteristics in units of lifetime
triplet = 0;
resolution = [30 10];
ring = [];
maxm = 10;

satv = [0 0.05*exp(log(1/0.05)*(0:9)/9)];


for jrad = 1:5
    for jcov = 1:1
        for jsat = 1:1

            over = radv(jrad);
            atf = [1.52 covv(jcov)];
            tic
            exc = DICExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focposexc, pow, atf, resolution, maxm);
            mdf = DICExc2MDF(exc, NA, n0, n, n1, focposdet, lamem, mag, av, zpin, atf, kappa, lt, pulse, satv(jsat), triplet);

            phi=0:pi/100:2*pi;

            FocusImage3D(exc.rho,exc.z-15,mdf.volx1+mdf.voly1);
            hold on
            FocusImage3D(exc.rho,exc.z-15,mdf.volx2+mdf.voly2);
            axis image
            hold off

            eval(['print -dpng -r300 DICFCSMDF0' mint2str(jrad,2) 'rad_' mint2str(jcov,2) 'cov_' mint2str(jsat,2) 'sat'])

        end
    end
end
