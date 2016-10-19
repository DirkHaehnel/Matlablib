if 1
    if 0
        clear fullres autotime radv covv jrad jcov jast w0 a0 dc vg v1 v2 mx1 my1 mx2 my2 wx1 wy1 wx2 wy2 amp1 amp2 w1 w2 a1 a2 rho z
        radv = [1.25 1.5 2 3 4]*1e3;
        covv = 0:10;
        astv = 0:0.1:0.3;
        jrad0=1, jcov0=1, jast0=1
    else
        load 2fFCSRadSatRef_full2
        jrad0=jrad, jcov0=jcov, jast0=jast
    end

    dist = 0.40;
    rhofield = [0 3];
    zfield = [10 20];
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
    resolution = [40 10];
    ring = [];
    maxm = 10;

    satv = [0 0.05*exp(log(1/0.05)*(0:9)/9)];
    %satv = 0;

    rmax = rhofield(end)/1.5;
    [x,y] = meshgrid(-rmax:0.05:rmax,-rmax:0.05:rmax);
    rr = sqrt(x.^2+y.^2);
    pp = angle(x+i*y);

    global pd
    pd = 1/5e-5;

    for jrad = 1:length(radv)
        for jcov = 1:length(covv)
            for jast = 1:length(astv)
                if jrad>jrad0 | (jrad==jrad0 & jcov>jcov0) | (jrad==jrad0 & jcov==jcov0 & jast>=jast0)
                    [jrad jcov jast]
                    over = [astv(jast) radv(jrad)];
                    atf = [1.52 covv(jcov)];
                    tic
                    for jsat = 1:length(satv)
                        [modres, autotime, exc, mdf] = DICFCS(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focposexc, ...
                            pow, lamem, mag, av, focposdet, zpin, atf, kappa, lt, pulse, satv(jsat), triplet, resolution, ring, maxm);
                        fullres(:,:,jrad,jcov,jast,jsat) = modres;
                        for jz=1:size(exc.z,2)
                            f1 = interp1(exc.rho(:,1),mdf.volx1(:,jz,1)+mdf.voly1(:,jz,1),rr,'cubic');
                            for j=1:maxm
                                f1 = f1 + interp1(exc.rho(:,1),mdf.volx1(:,jz,j+1)+mdf.voly1(:,jz,j+1),rr,'cubic').*cos(j*pp) + ...
                                    interp1(exc.rho(:,1),mdf.volx1(:,jz,maxm+1+j)+mdf.voly1(:,jz,maxm+1+j),rr,'cubic').*sin(j*pp);
                            end
                            f2 = interp1(exc.rho(:,1),mdf.volx2(:,jz,1)+mdf.voly2(:,jz,1),rr,'cubic');
                            for j=1:maxm
                                f2 = f2 + interp1(exc.rho(:,1),mdf.volx2(:,jz,j+1)+mdf.voly2(:,jz,j+1),rr,'cubic').*cos(j*pp) + ...
                                    interp1(exc.rho(:,1),mdf.volx2(:,jz,maxm+1+j)+mdf.voly2(:,jz,maxm+1+j),rr,'cubic').*sin(j*pp);
                            end

                            [mx1(jz,jrad,jcov,jast,jsat),my1(jz,jrad,jcov,jast,jsat),wx1(jz,jrad,jcov,jast,jsat),wy1(jz,jrad,jcov,jast,jsat),amp1(jz,jrad,jcov,jast,jsat)] = Gauss2D(x(1,:),y(:,1),f1+0.01*max(f1(:))*rand(size(f1)));
                            [mx2(jz,jrad,jcov,jast,jsat),my2(jz,jrad,jcov,jast,jsat),wx2(jz,jrad,jcov,jast,jsat),wy2(jz,jrad,jcov,jast,jsat),amp2(jz,jrad,jcov,jast,jsat)] = Gauss2D(x(1,:),y(:,1),f2+0.01*max(f2(:))*rand(size(f1)));
                        end

                        rho=exc.rho(:,1); z=exc.z(1,:);
                        w1(jrad,jcov,jast,jsat) = simplex('LaserBeamFit',0.4,[],[],[],[],lamex/n,z-15,sqrt((wx1(:,jrad,jcov,jast,jsat).^2+wy1(:,jrad,jcov,jast,jsat).^2)/2),1);
                        w2(jrad,jcov,jast,jsat) = simplex('LaserBeamFit',0.4,[],[],[],[],lamex/n,z-15,sqrt((wx2(:,jrad,jcov,jast,jsat).^2+wy2(:,jrad,jcov,jast,jsat).^2)/2),1);
                        a1(jrad,jcov,jast,jsat) = simplex('D1CEF',0.4,[],[],[],[],av/mag,lamem/n,z-15,amp1(:,jrad,jcov,jast,jsat));
                        a2(jrad,jcov,jast,jsat) = simplex('D1CEF',0.4,[],[],[],[],av/mag,lamem/n,z-15,amp2(:,jrad,jcov,jast,jsat));

                        para = simplex('GaussFcs',[0.4 0.1],[0 0],[],[],[],av/mag,[lamex lamem]/n,dist,...
                            autotime,fullres(:,1:2,1,jrad,jcov,jast,jsat),fullres(:,3,1,jrad,jcov,jast,jsat)/2);
                        w0(jrad,jcov,jast,jsat) = para(1);
                        a0(jrad,jcov,jast,jsat) = para(2);

                        dc(jrad,jcov,jast,jsat) = 1/pd/5e-5;

                        vg(jrad,jcov,jast,jsat) = GaussDetectionVolume(para(1), para(2), av/mag, [lamex lamem]/n);
                        v1(jrad,jcov,jast,jsat) = GaussDetectionVolume(w1(jrad,jcov,jast,jsat), a1(jrad,jcov,jast,jsat), av/mag, [lamex lamem]/n);
                        v2(jrad,jcov,jast,jsat) = GaussDetectionVolume(w2(jrad,jcov,jast,jsat), a2(jrad,jcov,jast,jsat), av/mag, [lamex lamem]/n);

                        save 2fFCSRadSatRef_full2 fullres autotime radv covv astv jrad jcov jast w0 a0 dc vg v1 v2 mx1 my1 mx2 my2 wx1 wy1 wx2 wy2 amp1 amp2 w1 w2 a1 a2 rho z

                    end
                    toc
                end
            end
        end
    end
end
