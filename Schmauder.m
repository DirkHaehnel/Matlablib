% program Ralf Schmauder

clear all
close all
rhofield = [0 1.5];
NA = 1.2;
dnp = 0.0186;
ng = 1.333;
% nb = 1.3336; 
% nb = 1.3618;
mag = 40;
over = [NA*164.5/mag*1e3 0 3e3];
focpos = 200;
zpin = 0;
atf = [];
kappa = 1;
lt = [];
sat = 0;
resolution = 30;
nbv = [ng 1.335:0.005:1.4]

if 0 % detection volume 1
    av = 70/2;
    lamex = 0.488;
    lamem = 0.525;
    zpin = 0; 
    phi0 = pi/2000:pi/1000:asin(NA/ng);    
        
    for jnb=1:length(nbv)
        nb = nbv(jnb)
        phi = asin(ng*sin(phi0)/nb);
        tmp = mean(tan(phi0)./tan(phi))*focpos;
        [volx, voly, rho, z, h0, h1, h2] = ...
            MDF(0, tmp + [-min([tmp,10]) 10], NA, ng, nb, nb, [], 0, [], lamex, over, focpos, lamem, mag, av, zpin, atf, kappa, lt, sat, 10);
        exc = volx(:,:,1)+voly(:,:,1); exc = exc/max(exc);
        z(exc<1e-3) = [];
        zfield = [z(1) z(end)]; 
        rhofield = [0 diff(zfield)/3];    
        [volx, voly, rho, z, h0, h1, h2] = ...
            MDF(rhofield, zfield, NA, ng, nb, nb, [], 0, [], lamex, over, focpos, lamem, mag, av, zpin, atf, kappa, lt, sat, resolution);
        [veffc(jnb), veff(jnb)] = DetectionVolume(rho,z,volx,voly);
        feld = [volx(:,:,1)+voly(:,:,1)];
        mpc(rho,z,feld,'lr');
        xlabel('\rho [\mum]');
        ylabel('\itz\rm [\mum]');
        title(['\itn_{sol}\rm = ' mnum2str(nb,1,3) ': \itV_{eff}\rm = ' mnum2str(veff(end),1,2) ' fl']);
        eval(['print -dpng -r100 Ralf488Refrac' mint2str(nb*1e3)])
    end
    save RalfVeff1 nbv veff veffc
        
end

if 0 % detection volume 2
    av = 90/2;
    lamex = 0.633;
    lamem = 0.670;
    zpin = 0; 
    phi0 = pi/2000:pi/1000:asin(NA/ng);    
    
    for jnb=1:length(nbv)
        nb = nbv(jnb)
        phi = asin(ng*sin(phi0)/nb);
        tmp = mean(tan(phi0)./tan(phi))*focpos;
        [volx, voly, rho, z, h0, h1, h2] = ...
            MDF(0, tmp + [-min([tmp,10]) 10], NA, ng, nb, nb, [], 0, [], lamex, over, focpos, lamem, mag, av, zpin, atf, kappa, lt, sat, 10);
        exc = volx(:,:,1)+voly(:,:,1); exc = exc/max(exc);
        z(exc<1e-3) = [];
        zfield = [z(1) z(end)]; 
        rhofield = [0 diff(zfield)/3];    
        [volx, voly, rho, z, h0, h1, h2] = ...
            MDF(rhofield, zfield, NA, ng, nb, nb, [], 0, [], lamex, over, focpos, lamem, mag, av, zpin, atf, kappa, lt, sat, resolution);
        [veffc(jnb), veff(jnb)] = DetectionVolume(rho,z,volx,voly);
        feld = [volx(:,:,1)+voly(:,:,1)];
        mpc(rho,z,feld,'lr');
        xlabel('\rho [\mum]');
        ylabel('\itz\rm [\mum]');
        title(['\itn_{sol}\rm = ' mnum2str(nb,1,3) ': \itV_{eff}\rm = ' mnum2str(veff(end),1,2) ' fl']);
        eval(['print -dpng -r100 Ralf633Refrac' mint2str(nb*1e3)])
    end
    save RalfVeff2 nbv veff veffc
    
end

if 0 % full calc. 1
    av = 70/2;
    lamex = 0.488;
    lamem = 0.525;
    zpin = 0; 
        
    for jnb=1:length(nbv)
        nb = nbv(jnb)
        [modres, autotime, veff, exc, rho, z, volx, voly] = ...
            NewFCS(NA, ng, nb, nb, [], 0, [], lamex, lamem, mag, focpos, av, over, zpin, atf, kappa, lt, sat, resolution);
        eval(['save Ralf488Refrac' mint2str(nb*1e3) ' modres autotime veff exc rho z volx voly'])
    end
end

if 0 % full calc. 2
    av = 90/2;
    lamex = 0.633;
    lamem = 0.670;
    zpin = 0; 
    
    for jnb=1:length(nbv)
        nb = nbv(jnb)
        [modres, autotime, veff, exc, rho, z, volx, voly] = ...
            NewFCS(NA, ng, nb, nb, [], 0, [], lamex, lamem, mag, focpos, av, over, zpin, atf, kappa, lt, sat, resolution);
        eval(['save Ralf633Refrac' mint2str(nb*1e3) ' modres autotime veff exc rho z volx voly'])
    end
end

if 0
    for jnb=1:length(nbv)
        nb = nbv(jnb)
        eval(['load Ralf488Refrac' mint2str(nb*1e3) ' modres autotime rho volx voly'])        
        auto1(:,jnb) = modres;
        vx1(jnb) = sum(sum(rho.*volx(:,:,1)));
        vy1(jnb) = sum(sum(rho.*voly(:,:,1)));        
        eval(['load Ralf633Refrac' mint2str(nb*1e3) ' modres autotime rho volx voly'])        
        auto2(:,jnb) = modres;
        vx2(jnb) = sum(sum(rho.*volx(:,:,1)));
        vy2(jnb) = sum(sum(rho.*voly(:,:,1)));        
    end
    nn=nbv(1):0.001:1.4;
    v1 = interp1([2*1.333-nbv(end:-1:2) nbv],[vx1(end:-1:2)+vy1(end:-1:2) vx1+vy1],nn,'cubic');
    v2 = interp1([2*1.333-nbv(end:-1:2) nbv],[vx2(end:-1:2)+vy2(end:-1:2) vx2+vy2],nn,'cubic');    
    save RalfResults nn v1 v2 autotime auto1 auto2 vx1 vx2 vy1 vy2 nbv
    
    return
    
    plot(nn,v1/max(v1),nn,v2/max(v2))
    axis([1.333 1.4 0.3 1])
    xlabel('refractive index'); ylabel('fluorescence intensity');
    legend({'488/525 nm & 70 \mum','633/670 nm & 90 \mum'})
    
    fout=fopen('RalfACF1.txt','w')
    for j=1:length(autotime) fprintf(fout,'%f\t', autotime(j)); for k=1:size(auto1,2) fprintf(fout,'%f\t', auto1(j,k)); end; fprintf(fout,'\n'); end
    fclose(fout)
    fout=fopen('RalfACF2.txt','w')
    for j=1:length(autotime) fprintf(fout,'%f\t', autotime(j)); for k=1:size(auto1,2) fprintf(fout,'%f\t', auto2(j,k)); end; fprintf(fout,'\n'); end
    fclose(fout)
end

if 1 % full calc. for measurements in RalfData
    load RalfData.mat
    load Ralf488ModelAstSat3mm
    nbv = nv;
    av = 70/2;
    lamex = 0.488;
    lamem = 0.525;
    astv = 0:0.1:0.5;
    atfv = -5:5;
    satv = 0:0.1:1;
    for jnb=5:length(nbv)
        nb = nbv(jnb);
        for jatf = 1:length(atfv)
            atf = [1.52 atfv(jatf)];
            for jast = 1:length(astv)
                over = [NA*164.5/mag*1e3 astv(jast) 10e3 1];
                for jsat = 1:length(satv)
                    sat = satv(jsat);
                    [modres(:,jatf,jast,jsat,jnb), autotime] = ...
                        NewFCS(NA, ng, nb, nb, [], 0, [], lamex, lamem, mag, focpos, av, over, zpin, atf, kappa, lt, sat, resolution);
                    eval(['save Ralf488ModelAstSat3mm nbv satv astv atfv modres autotime'])
                end
            end
        end
    end
end

if 1 % fitting of experiment
    load Ralf488ModelAstSat3mm
    autotime0 = autotime;
    load RalfData.mat
    for j=1:length(nv)
        for jatf=1:length(atfv)
            for jast=1:length(astv)
                for jsat=1:length(satv)
                    p(:,jatf,jast,jsat,j) = simplex('affineexpfit',[1/5 1e1],[0 0],[],[],[], autotime, auto(:,j), autotime0, modres(:,jatf,jast,jsat,j), 1);
                    err(jatf,jast,jsat,j) = affineexpfit(p(:,jatf,jast,jsat,j), autotime, auto(:,j), autotime0, modres(:,jatf,jast,jsat,j));
                end
            end
        end
    end
    save RalfFitResultsAstSat3mm autotime auto autotime0 modres nbv astv atfv satv p err
end

