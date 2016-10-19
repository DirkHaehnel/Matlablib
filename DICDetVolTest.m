if 1
    if 0
        clear fullres autotime radv covv jrad jcov jast w0 a0 dc vg v1 v2 mx1 my1 mx2 my2 wx1 wy1 wx2 wy2 amp1 amp2 w1 w2 a1 a2 rho z
        radv = [1.25 1.5 2 3 4]*1e3;
        covv = 0:2.5:10;
        jrad0=1, jcov0=1,
    else
        load DICDetVolRes
        jrad0=jrad; jcov0=jcov;
    end

    dist = 0.40;
    rhofield = [0 8];
    zfield = [5 25];
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
    sat = 0;

    global pd
    pd = 1/5e-5;

    for jrad = 1:length(radv)
        for jcov = 1:length(covv)
            if jrad>jrad0 || (jrad==jrad0 && jcov>=jcov0)
                [jrad jcov]
                over = [0 radv(jrad)];
                atf = [1.52 covv(jcov)];
                tic
                [modres, autotime, exc, mdf] = DICFCS(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focposexc, ...
                    pow, lamem, mag, av, focposdet, zpin, atf, kappa, lt, pulse, sat, triplet, resolution, ring, maxm);
                fullres(:,:,jrad,jcov) = modres;

                para = simplex('GaussFcs',[0.4 0.1],[0 0],[],[],[],av/mag,[lamex lamem]/n,dist,...
                    autotime,fullres(:,1:2,jrad,jcov),fullres(:,3,jrad,jcov)/2);
                w0(jrad,jcov) = para(1);
                a0(jrad,jcov) = para(2);

                vg(jrad,jcov) = GaussDetectionVolume(para(1), para(2), av/mag, [lamex lamem]/n);
                tst=[]; 
                for j=5:10 
                    ind = abs(exc.z(1,:)-15)<=j; 
                    tst(j) = DetectionVolume(exc.rho(:,ind), exc.z(:,ind), mdf.volx1(:,ind,:), mdf.voly1(:,ind,:)); 
                end
                para = simplex('ExpFun',1,0,[],[],[],5:10,15-tst(5:10));
                [err,c] = ExpFun(para,5:10,15-tst(5:10));
                v1(jrad,jcov) = 15-c(1);
                tst=[]; 
                for j=5:10 
                    ind = abs(exc.z(1,:)-15)<=j; 
                    tst(j) = DetectionVolume(exc.rho(:,ind), exc.z(:,ind), mdf.volx2(:,ind,:), mdf.voly2(:,ind,:)); 
                end
                para = simplex('ExpFun',1,0,[],[],[],5:10,15-tst(5:10));
                [err,c] = ExpFun(para,5:10,15-tst(5:10));
                v2(jrad,jcov) = 15-c(1);
                
                toc
                
                save DICDetVolRes fullres autotime radv covv jrad jcov w0 a0 vg v1 v2 jrad jcov
            end
        end
    end
end
