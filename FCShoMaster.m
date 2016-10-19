if 1
    if 1
        clear all
        radv = [1.25 1.5 2 3 4]*1e3;
        covv = 0:10;
        astv = 0:0.1:0.3;
        jrad0=1, jcov0=1, jast0=1
    else
        load FCShoRadSatRef
        jrad0=jrad, jcov0=jcov, jast0=jast
    end

    dist = 0.40;
    rhofield = [0 4];
    zfield = [5 25];
    NA = 1.2;
    fd = 3e3;
    n0 = 1.333;
    n = n0;
    n1 = n0;
    d0 = [];
    d = 0;
    d1 = [];
    lamex = 0.64;
    focpos = 15;
    pow = 1;
    lamem = 0.67;
    mag = 60;
    av = 100;
    zpin = 0e3;
    kappa = 0;
    lt = [];
    pulse = [0.05 25]/2; % laser characteristics in units of lifetime
    triplet = 0;
    resolution = [40 10];
    ring = [];
    maxm = 10;
    ord = 4;

    satv = [0 0.05*exp(log(1/0.05)*(0:9)/9)];
    %satv = 0;

    rmax = rhofield(end)/1.5;
    [x,y] = meshgrid(-rmax:0.05:rmax,-rmax:0.05:rmax);
    rr = sqrt(x.^2+y.^2);
    pp = angle(x+i*y);

    for jrad = 1:length(radv)
        for jcov = 1:length(covv)
            for jast = 1:length(astv)
                if jrad>jrad0 || (jrad==jrad0 && jcov>jcov0) || (jrad==jrad0 && jcov==jcov0 && jast>=jast0)
                    [jrad jcov jast]
                    over = [astv(jast) radv(jrad)];
                    atf = [1.52 covv(jcov)];
                    tic
                    exc = GaussExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, atf, resolution, ring, maxm);
                    for jsat = 1:length(satv)
                        sat = satv(jsat);
                        mdf = GaussExc2MDF(exc, NA, n0, n, n1, focpos, lamem, mag, av, zpin, atf, kappa, lt, pulse, sat, triplet);
                        [veff(:,jrad,jcov,jast,jsat), intens(:,jrad,jcov,jast,jsat)] = DetectionVolume(exc.rho,exc.z,mdf.volx+mdf.voly,repmat(mdf.volx+mdf.voly,[1 1 1 3]));
                        [modres, autotime] = FCSho(exc.rho, exc.z, mdf.volx+mdf.voly, 0, 4);
                        fullres(:,:,jrad,jcov,jast,jsat) = modres;
                        for jz=1:size(exc.z,2)
                            f = interp1(exc.rho(:,1),mdf.volx(:,jz,1)+mdf.voly(:,jz,1),rr,'cubic');
                            for j=1:maxm
                                f = f + interp1(exc.rho(:,1),mdf.volx(:,jz,j+1)+mdf.voly(:,jz,j+1),rr,'cubic').*cos(j*pp) + ...
                                    interp1(exc.rho(:,1),mdf.volx(:,jz,maxm+1+j)+mdf.voly(:,jz,maxm+1+j),rr,'cubic').*sin(j*pp);
                            end
                            [mx(jz,jrad,jcov,jast,jsat),my(jz,jrad,jcov,jast,jsat),wx(jz,jrad,jcov,jast,jsat),wy(jz,jrad,jcov,jast,jsat),amp(jz,jrad,jcov,jast,jsat)] = Gauss2D(x(1,:),y(:,1),f+0.01*max(f(:))*rand(size(f)));
                        end
                        rho=exc.rho(:,1); z=exc.z(1,:);
                        w(jrad,jcov,jast,jsat) = simplex('LaserBeamFit',0.4,[],[],[],[],lamex/n,z-15,sqrt(wx(:,jrad,jcov,jast,jsat).*wy(:,jrad,jcov,jast,jsat)),1);
                        a(jrad,jcov,jast,jsat) = simplex('D1CEF',0.4,[],[],[],[],av/mag,lamem/n,z-15,amp(:,jrad,jcov,jast,jsat));

                        para = simplex('Riglerho',[100 1e3],[0 0],[],[],[],autotime,modres(:,1),2);
                        a2(jrad,jcov,jast,jsat) = para(1);
                        b2(jrad,jcov,jast,jsat) = para(2);
                        para = simplex('Riglerho',para,[0 0],[],[],[],autotime,modres(:,2),3);
                        a3(jrad,jcov,jast,jsat) = para(1);
                        b3(jrad,jcov,jast,jsat) = para(2);
                        para = simplex('Riglerho',para,[0 0],[],[],[],autotime,modres(:,3),4);                        
                        a4(jrad,jcov,jast,jsat) = para(1);
                        b4(jrad,jcov,jast,jsat) = para(2);
                        
                        save FCShoRadSatRef 
                    end
                    toc
                end
            end
        end
    end
end
