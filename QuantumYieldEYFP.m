% Program for Frank's cavity project

% clear all

if 1
    load metals
    lamem = 0.535;
    lamex = 0.47;

    nglass = 1.52; % glass
    nsio = 1.46;
    npva = 1.52;
    nkleber = 1.56;
    nagex = silver(wavelength==1e3*lamex);
    nagem = silver(wavelength==1e3*lamem);
    dagbottom = 0.03;
    dagtop = 0.06;
    dsio = 0.02;
    dpva = 0.05;
    dkleberv = 0.1:0.001:1;

    rhofield = [0 1];
    NA = 1.46;
    fd = 164.5e3/100;
    over = 3e3;
    focpos = 0;
    mag = 100;
    av = 75;
    zpin = 0;
    atf = [];
    kappa = 1;
    for jd = 168:length(dkleberv)
        dkleber = dkleberv(jd);
        resolution = [20 round(lamex/0.01)]; % 10 nm along z

        n0 = [nglass nagex nsio]; d0 = [dagbottom dsio];
        n = npva; d = dpva;
        n1 = [nkleber nagex nglass]; d1 = [dkleber dagtop];
        
        exc = GaussExc(rhofield, [0, dpva], NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, atf, resolution);
        
        n0 = [nglass nagem nsio];
        n1 = [nkleber nagem nglass];
        
        mdf = GaussExc2MDF(exc, NA, n0, n, n1, focpos, lamem, mag, av, zpin, atf, kappa, 1);
        tau = GaussExc2MDF(exc, NA, n0, n, n1, focpos, lamem, mag, av, zpin, atf, kappa, 2);
        
        res(jd) = 4/3*n*squeeze(sum(sum(exc.rho.*(tau.volx(:,:,1)+tau.voly(:,:,1)))))/squeeze(sum(sum(exc.rho.*(mdf.volx(:,:,1)+mdf.voly(:,:,1)))));
        brightness(jd) = squeeze(sum(sum(exc.rho.*(mdf.volx(:,:,1)+mdf.voly(:,:,1)))))*diff(exc.z(1,1:2));
        jd
    end
    save QY_EYFP
end

