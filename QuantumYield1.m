% Program for the quantum yield project

% clear all

if 1
    load metals
    lamex = 0.635;
    lamem = 0.685;

    ng = 1.52;
    n = 1.4327;
    nauex = gold(wavelength==1e3*lamex);
    nauem = gold(wavelength==1e3*lamem);
    d0v = [0.01:0.02:0.07];
    d1 = 0.1;
    rhofield = [0 1];
    NA = 1.4;
    fd = 3e3;
    dv = 0.0025:0.005:0.5;
    over = 3e3;
    focposv = [-1 0 1];
    mag = 60;
    av = 50;
    zpin = 0;
    atf = [];
    kappa = 1;
    resolution = 50;
    delta = 1e-3;
    for jfp = 1:length(focposv)
        focpos = focposv(jfp)
        for jd0 = 1:length(d0v)
            d0 = d0v(jd0);
            for j=1:length(dv)
                resolution = [20 round(lamex/(dv(j)/100))];
                exc = GaussExc(rhofield, [0, dv(j)], NA, fd, [ng nauex], n, [nauex ng], d0, dv(j), d1, lamex, over, focpos, atf, resolution);
                mdf = GaussExc2MDF(exc, NA, [ng nauem], n, [nauem ng], focpos, lamem, mag, av, zpin, atf, kappa, 1);
                tau = GaussExc2MDF(exc, NA, [ng nauem], n, [nauem ng], focpos, lamem, mag, av, zpin, atf, kappa, 2);
                % eval(['save QYGlycerol_' mint2str(jfp,1) '_' mint2str(jd0,1) '_' mint2str(j,3) ' exc mdf tau']);
                res(j,jd0,jfp) = 4/3*n*squeeze(sum(sum(exc.rho.*(tau.volx(:,:,1)+tau.voly(:,:,1)))))/squeeze(sum(sum(exc.rho.*(mdf.volx(:,:,1)+mdf.voly(:,:,1)))));
                brightness(j,jd0,jfp) = squeeze(sum(sum(exc.rho.*(mdf.volx(:,:,1)+mdf.voly(:,:,1)))))*diff(exc.z(1,1:2));
                j
            end
            eval(['save QYGlycerol dv brightness res']);
        end
    end

    return

    di=0.0025:0.0001:dv(end);
    plot(di,interp1(dv,brightness(:,:,1),di,'spline'))
    plot(di,interp1(dv,res(:,:,1),di,'spline'))

    ind = dv>0.15 & dv<0.25;
    plot([0 d0v],[1 min(res(ind,:,1))])
end

if 0
    load metals
    lamex = 0.635;
    lamem = 0.685;

    ng = 1.52;
    n = 1.4327;
    nauex = gold(wavelength==1e3*lamex);
    nauem = gold(wavelength==1e3*lamem);
    d0v = [0.01:0.02:0.07];
    d1 = 0.1;
    rhofield = [0 1];
    NA = 1.4;
    fd = 3e3;
    dv = 0.0025:0.005:0.5;
    over = 3e3;
    focposv = [-1 0 1];
    mag = 60;
    av = 50;
    zpin = 0;
    atf = [];
    kappa = 1;
    resolution = 50;
    delta = 1e-3;
    for jfp = 1:length(focposv)
        focpos = focposv(jfp)
        for jd0 = 1:length(d0v)
            d0 = d0v(jd0);
            for j=1:length(dv)
                reflex(j,jd0,jfp) = QuantumYieldReflex(1.4, 3e3, 60, [ng nauex n nauex ng], [d0, dv(j), d1], lamex, over, focpos);
                j
            end
            eval(['save QYGlycerolReflex dv reflex']);
        end
    end
end