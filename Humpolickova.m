if 1
    rhofield = [0 3];
    zfield = [10 20];
    NA = 1.2;
    fd = 3e3;
    n0 = 1.333;
    n = n0;
    n1 = n0;
    d0 = [];
    d = 0;
    d1 = [];
    lamex = 0.633;
    focpos = 15;
    lamem = 0.67;
    mag = 60;
    av = 45;
    zpin = 0e3;
    kappa = 0;
    lt = [];
    pulse = []; % laser characteristics in units of lifetime
    triplet = 0;
    resolution = [30 10];
    ring = [];
    maxm = 2;

    rmax = rhofield(end)/1.5;
    [x,y] = meshgrid(-rmax:0.05:rmax,-rmax:0.05:rmax);
    rr = sqrt(x.^2+y.^2);
    pp = angle(x+i*y);

    over = 2e3;
    atf = [];
    exc = GaussExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, atf, resolution);
    mdf = GaussExc2MDF(exc, NA, n0, n, n1, focpos, lamem, mag, av, zpin);

    for jz=1:size(exc.z,2)
        ff = interp1(exc.rho(:,1),mdf.volx(:,jz,1)+mdf.voly(:,jz,1),rr,'cubic');
        for j=1:maxm
            ff = ff + interp1(exc.rho(:,1),mdf.volx(:,jz,j+1)+mdf.voly(:,jz,j+1),rr,'cubic').*cos(j*pp) + ...
                interp1(exc.rho(:,1),mdf.volx(:,jz,maxm+1+j)+mdf.voly(:,jz,maxm+1+j),rr,'cubic').*sin(j*pp);
        end
        [mx(jz),my(jz),wx(jz),wy(jz),amp(jz)] = Gauss2D(x(1,:),y(:,1),ff+0.01*max(ff(:))*rand(size(ff)));
    end

    rho=exc.rho(:,1); z=exc.z(1,:);
    ww = simplex('LaserBeamFit',0.4,[],[],[],[],lamex/n,z(abs(z-15)<=1.2)-15,sqrt((wx(abs(z-15)<=1.2).^2+wy(abs(z-15)<=1.2).^2)/2),[],1);
    aa = simplex('D1CEF',0.4,[],[],[],[],av/mag,lamem/n,z(abs(z-15)<=1.2)-15,amp(abs(z-15)<=1.2),1);

end
