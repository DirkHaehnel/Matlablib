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
    d0 = 0.05;
    d = 0.1;
    d1 = 0.1;
    rhofield = [0 0.4];
    NA = 1.4;
    fd = 3e3;
    dv = 0.0025:0.005:0.5;
    over = 3e3;
    focpos = 0;
    mag = 60;
    av = 50;
    zpin = 0;
    atf = [];
    kappa = 1;
    resolution = 50;
    delta = 1e-3;
	resolution = [20 round(lamex/(d/100))];
    exc = GaussExc(rhofield, [0, d], NA, fd, [ng nauex], n, [nauex ng], d0, d, d1, lamex, over, focpos, atf, resolution);
    mdf = GaussExc2MDF(exc, NA, [ng nauem], n, [nauem ng], focpos, lamem, mag, av, zpin, atf, kappa, 1);
    tau = GaussExc2MDF(exc, NA, [ng nauem], n, [nauem ng], focpos, lamem, mag, av, zpin, atf, kappa, 2);

    subplot(3,1,1);
    FocusImage2D(exc.rho,exc.z,cat(3,exc.fxc,exc.fxs));
    colorbar('h');
    subplot(3,1,2);
    FocusImage2D(exc.rho,exc.z,mdf.volx+mdf.voly);
    colorbar('h');
    subplot(3,1,3);
    FocusImage2D(exc.rho,exc.z,ones(size(exc.rho,1),1)*(1./tau.qv));
    colorbar('h')
    
    res = 4/3*n*squeeze(sum(sum(exc.rho.*(tau.volx(:,:,1)+tau.voly(:,:,1)))))/squeeze(sum(sum(exc.rho.*(mdf.volx(:,:,1)+mdf.voly(:,:,1)))));
	brightness = squeeze(sum(sum(exc.rho.*(mdf.volx(:,:,1)+mdf.voly(:,:,1)))))*diff(exc.z(1,1:2));

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