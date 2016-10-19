clear all
close all

transmission = [821.36 801.77 782.72 764.22 746.27 731.52 711.48 694.91 680.07 664.05 646.91 633.79 622.61 611.20 598.70 585.89 573.48 562.28 552.75 543.02 534.34 525.05 514.84 509.14 502.13 495.67 489.76 484.39 479.57 475.30 471.57 468.39 465.76 463.68 462.14 461.15 460.71 460.82 461.47];
nglass = 1.52; % glass
nsio = 1.46;
npva = 1.52;
nkleber = 1.56;
dagbottom = 0.03;
dagtop = 0.06;
dsio = 0.02;
dpva = 0.05;
dkleberv = 0.025:0.001:0.15;

if 0 % Lifetime fitting
    load c:\Joerg\Doc\Nanophotonics&Optics\FrankSchleifenbaum\PumpProbe\EYFP\irf.mat % irf from 6.11.07

    cd c:\Joerg\Doc\Nanophotonics&Optics\FrankSchleifenbaum\PumpProbe\EYFP
    fnames = dir('*.txt');
    for j=1:length(fnames)
        if strfind(fnames(j).name,'_IRF')
            tmp = j;
        end
    end
    fnames(tmp) = [];
    for j=1:length(fnames)
        fin = fopen(fnames(j).name,'r');
        for k=1:3
            fgetl(fin);
        end
        tmp = fscanf(fin,'%f %f\n',inf);
        fclose(fin);
        y(:,:,j) = [tmp(1:2:end) tmp(2:2:end)];
    end

    ind = irf(:,1)>35 & irf(:,1)<60;
    dt = mean(diff(irf(:,1))); % TCSPC channel
    p = 100/dt; % periodicity of TCSPC excitation
    tp = (1:p)';
    tau = [0.5 1 2];
    irf = circshift(irf(ind,2),-50);
    for j=1:size(y,3)
        close all

        [c, offset, A, tau] = fluofit(irf, y(ind,2,j), p, dt, tau);
        res(j).c = c;
        res(j).offset = offset;
        res(j).A = A;
        res(j).tau = tau;

        eval(['print -dpng -r300 ' fnames(j).name(1:end-4)])
    end
    for j=1:39 tauav(j)=res(j).A'*res(j).tau; end
    save EYFPData res fnames tauav
end

if 1 % Cavity lifetime modelling
    load SilverDmitri
    lamem = 0.535;
    lamex = 0.47;

    nagex = interp1(SilverDmitri(:,1),SilverDmitri(:,2),1e3*lamex);
    nagem = interp1(SilverDmitri(:,1),SilverDmitri(:,2),1e3*lamem);

    rhofield = [0 1];
    NA = 1.3;
    fd = 164.5e3/100;
    over = 3e3;
    focpos = 0;
    mag = 300;
    av = 75;
    zpin = 0;
    atf = [];
    kappa = 1;
    resolution = [20 floor(lamex/0.005)]; % 5 nm resolution
    for jd = 1:length(dkleberv)
        dkleber = dkleberv(jd);

        n0 = [nglass nagex nsio]; d0 = [dagbottom dsio];
        n = npva; d = dpva;
        n1 = [nkleber nagex nglass]; d1 = [dkleber dagtop];
        
        exc = GaussExc(rhofield, [0, dpva], NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, atf, resolution);
        
        n0 = [nglass nagem nsio];
        n1 = [nkleber nagem nglass];
        
        mdf = GaussExc2MDF(exc, NA, n0, n, n1, focpos, lamem, mag, av, zpin, atf, kappa, 1);
        tau = GaussExc2MDF(exc, NA, n0, n, n1, focpos, lamem, mag, av, zpin, atf, kappa, 2);
        
        lifetime0(jd) = 4/3*n*squeeze(sum(sum(exc.rho.*(tau.volx(:,:,1)+tau.voly(:,:,1)))))/squeeze(sum(sum(exc.rho.*(mdf.volx(:,:,1)+mdf.voly(:,:,1)))));
        brightness(jd) = squeeze(sum(sum(exc.rho.*(mdf.volx(:,:,1)+mdf.voly(:,:,1)))))*diff(exc.z(1,1:2));
        jd
    end
    save QY_EYFP lifetime0 brightness
end

if 1 % Cavity transmission
    load SilverDmitri
    lamexv = 450:850;

    clear tp rp ts rs
    for jd = 1:length(dkleberv)
        dkleber = dkleberv(jd);
        for jl = 1:length(lamexv)
            nagex = interp1(SilverDmitri(:,1),SilverDmitri(:,2),lamexv(jl));
            n0 = [nglass nagex nsio]; d0 = [dagbottom dsio];
            n = npva; d = dpva;
            n1 = [nkleber nagex nglass]; d1 = [dkleber dagtop];
            [rp(jl,jd), rs(jl,jd), tp(jl,jd), ts(jl,jd)] = Fresnel(nglass,[n0 n n1],1e3*[d0 d d1]/lamexv(jl)*2*pi);
        end
        jd
    end
    save TransmissionModel dkleberv lamexv rp rs tp ts
    pcolor((dsio+dpva+dkleberv)*1e3,lamexv,abs(tp).^2); shading interp    
    xlabel('cavity thickness (nm)'); ylabel('wavelength (nm)');
end

if 1 % data fitting
    load EYFPData
    load QY_EYFP
    load TransmissionModel

    for jd=1:size(tp,2)
        maxlam(jd) = lamexv(abs(tp(:,jd))==max(abs(tp(:,jd))));
    end
    ct = polyfit(maxlam,1e3*(dsio+dpva+dkleberv),1);
    thickness = polyval(ct,transmission);

    for k=1:50
        delta = (k-1)/10;
        p = Simplex('AffineFit',1,1,1,[],[],thickness', 1./tauav', 1e3*(dsio+dpva+dkleberv)' + delta, 1./lifetime0',[],1);
        [err(k), c, z] = AffineFit(p,thickness', 1./tauav', 1e3*(dsio+dpva+dkleberv)' + delta, 1./lifetime0',[],1);
    end
    [k,k] = min(err);
    delta = (k-1)/10;
    p = Simplex('AffineFit',1,1,1,[],[],thickness', 1./tauav', 1e3*(dsio+dpva+dkleberv)' + delta, 1./lifetime0',[],1);
    [err(k), c, z] = AffineFit(p,thickness', 1./tauav', 1e3*(dsio+dpva+dkleberv)' + delta, 1./lifetime0',[],1);
    plot(thickness,1./z,thickness,tauav,'o')
    xlabel('cavity thickness (nm)'); ylabel('average lifetime (ns)');
    
    quantumyield_EYFP = c(2)/sum(c);
    lifetime_EYFPwater = nkleber/1.33/sum(c);
end

