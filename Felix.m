if 1
    NA = 1.4;
    n0 = 1.51;
    n = 1.33;
    n1 = n;
    d0 = [];
    d = [];
    d1 = [];
    lamex = 0.635;
    lamem = 0.670;
    mag = 100*175/180;
    av = 100/2;
    kappa = 0;
    atf = [];
    lt = [];
    
    load felixexp
    dd = 209e3;
    p1 = mean(p1(3:4)); p2 = mean(p2(3:4));
    p1 = 2*9.3*p1; p2 = 2*9.3*p2; % one camera pixel = 9.3 mum
    w0 = simplex('BeamWaistFun',p1,0,min([p1,p2]),[],500,p1,p2,dd,lamex);
    z0 = sqrt((p1/w0)^2-1)/lamex*pi*w0^2;
    over = [1.8e3*NA z0 w0]
    
    zpv = -30e3:1e3:15e3;
    fpv = 5:25;
    for j = 1:length(zpv)
        for k=1:length(fpv)
            zpin = zpv(j);
            focpos = fpv(k);
            [zpin focpos]
            [modres(:,:,j,k), autotime, veff(:,j,k), exc(:,j,k)] = ...
                NewFCS(NA, n0, n, n1, d0, d, d1, lamex, lamem, mag, focpos, av, over, zpin, kappa, atf, lt);
        end
    end
    
%    save felixmod NA n0 n n1 d0 d d1 lamex lamem mag av over kappa atf lt modres autotime veff exc zpv fpv
end

if 1
    res = res3;
    t = res.autotime(10:end-10);
    x = mean(res.auto(10:end-10,:)')';
    tmp = []; 
    for j=1:size(modres,3)
        for k=1:size(modres,4)
            p(j,k) = simplex('affinefit',5e-5,0,[],[],[],t,x,autotime,modres(:,2,j,k)); 
            [err(j,k) tmp, z(:,j,k)] = affinefit(p(j,k),t,x,autotime,modres(:,2,j,k)); 
        end
    end
    diffusion = 5e-5./p*1e-8
end

if 0
    tm = 28.8;
    rad = 17.971e-9;
    natconst
    temp = 0:10:100;
    visc = [1.792 1.308 1.005 0.8007 0.6560 0.5494 0.4688 0.4061 0.3565 0.3165 0.2838]*1e-3;
    eta = interp1(temp,visc,tm,'cubic');
    tm = tm + 273.15;
    plot(zpv,diffusion,zpv,ones(size(zpv))*BoltzmannConstant*tm/6/pi/eta/rad*1e4)
end
