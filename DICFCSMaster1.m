if 1
    if 1
        clear fullres autotime radv covv jrad jcov jtrp w0 a0 dc vg v1 v2 mx1 my1 mx2 my2 wx1 wy1 wx2 wy2 amp1 amp2 w1 w2 a1 a2 rho z
        radv = [1.25 1.5 2 3 4]*1e3;
        covv = 0; %0:5:10;
        trpv = [0 0.5 1 2 4 8];
        jrad0=1, jcov0=1, jtrp0=1
    else
        load 2fFCSRadSatCovTrp
        jrad0=jrad, jcov0=jcov, jtrp0=jtrp
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
            for jtrp = 1:length(trpv)
                if jrad>jrad0 | (jrad==jrad0 & jcov>jcov0) | (jrad==jrad0 & jcov==jcov0 & jtrp>=jtrp0)
                    [jrad jcov jtrp]
                    over = radv(jrad);
                    atf = [1.52 covv(jcov)];
                    tic
                    for jsat = 1:length(satv)
                        triplet = trpv(jtrp);
                        [modres, autotime, exc, mdf] = DICFCS(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focposexc, ...
                            pow, lamem, mag, av, focposdet, zpin, atf, kappa, lt, pulse, satv(jsat), triplet, resolution, ring, maxm);
                        fullres(:,:,jrad,jcov,jtrp,jsat) = modres;
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

                            [mx1(jz,jrad,jcov,jtrp,jsat),my1(jz,jrad,jcov,jtrp,jsat),wx1(jz,jrad,jcov,jtrp,jsat),wy1(jz,jrad,jcov,jtrp,jsat),amp1(jz,jrad,jcov,jtrp,jsat)] = Gauss2D(x(1,:),y(:,1),f1+0.01*max(f1(:))*rand(size(f1)));
                            [mx2(jz,jrad,jcov,jtrp,jsat),my2(jz,jrad,jcov,jtrp,jsat),wx2(jz,jrad,jcov,jtrp,jsat),wy2(jz,jrad,jcov,jtrp,jsat),amp2(jz,jrad,jcov,jtrp,jsat)] = Gauss2D(x(1,:),y(:,1),f2+0.01*max(f2(:))*rand(size(f1)));
                        end

                        rho=exc.rho(:,1); z=exc.z(1,:);

                        pd = 1/5e-5;
                        para = simplex('GaussFcs',[0.4 0.1],[0 0],[],[],[],av/mag,[lamex lamem]/n,dist,...
                            autotime,fullres(:,1:2,jrad,jcov,jtrp,jsat),fullres(:,3,jrad,jcov,jtrp,jsat)/2);
                        w0(jrad,jcov,jtrp,jsat) = para(1);
                        a0(jrad,jcov,jtrp,jsat) = para(2);

                        dc(jrad,jcov,jtrp,jsat) = 1/pd/5e-5;

                        save 2fFCSRadSatCovTrp fullres autotime radv covv trpv jrad jcov jtrp w0 a0 dc mx1 my1 mx2 my2 wx1 wy1 wx2 wy2 amp1 amp2 rho z

                    end
                    toc
                end
            end
        end
    end
end
