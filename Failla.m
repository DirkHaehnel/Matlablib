% clear all

if 1
    load metals
    lamex = 0.488;
    lamem = 0.532;

	nw = 1.33;
    ng = 1.52;
    n = 1.56;
    nauex = silver(wavelength==1e3*lamex);
    nauem = silver(wavelength==1e3*lamem);
    d0 = 0.03;
    dv = 0.08:0.001:0.15;
    d1 = 0.06;
    rhofield = [0 0.4];
    NA = 1.2;
    fd = 1.65e3;
    over = 3e3;
    focpos = 0;
    mag = 100;
    av = 50;
    zpin = 0;
    atf = [];
    kappa = 1;
    delta = 1e-3;

    res = zeros(size(dv));
    brightness = res;
    for j=1:length(dv)
        d = dv(j);
        resolution = [20 round(lamex/(d/100))];
        exc = GaussExc(rhofield, [0 d], NA, fd, [ng nauex], n, [nauex ng], d0, d, d1, lamex, over, focpos, atf, resolution);
        mdf = GaussExc2MDF(exc, NA, [ng nauem], n, [nauem ng], focpos, lamem, mag, av, zpin, atf, kappa, 1);
        tau = GaussExc2MDF(exc, NA, [ng nauem], n, [nauem ng], focpos, lamem, mag, av, zpin, atf, kappa, 2);

        res(j) = 4/3*n*squeeze(sum(sum(exc.rho.*(tau.volx(:,:,1)+tau.voly(:,:,1)))))/squeeze(sum(sum(exc.rho.*(mdf.volx(:,:,1)+mdf.voly(:,:,1)))));
        brightness(j) = squeeze(sum(sum(exc.rho.*(mdf.volx(:,:,1)+mdf.voly(:,:,1)))))*diff(exc.z(1,1:2));
        j
    end
    save Failla
    plot(dv*1e3,1./res); axis([80 150 0 4]); xlabel('cavity height (nm)'); ylabel('emission rate (\itk\rm_0)');
end


return

if mod(1e3*d,10)==0
        FocusImage2D(exc.rho*1e3,exc.z*1e3,cat(3,exc.fxc,exc.fxs));
    eval(['print -dpng -r300 ''FrankCavityIntensity' mint2str(d*1e3,3) ' nm'''])
end


return

for j=1:size(exc.z,2)
    [lvd,lvu,lpd,lpu,qvd,qvu,qpd,qpu,qv,qp] = LifetimeL(exc.z(1,j)/lamem*2*pi, [ng nauex], n, [nauex ng], d0/lamem*2*pi, d/lamem*2*pi, d1/lamem*2*pi);
    taup(j) = 4/3*n/(qpu+qpd);
    tauv(j) = 4/3*n/(qvu+qvd);
end
q=1-1e-2; % quantum yield of photochemsitry
plot(exc.z(1,:)*1e3,q./(q+(1-q)./taup),exc.z(1,:)*1e3,q./(q+(1-q)./tauv));
xlabel('cavity height (nm)'); ylabel('quantum yield');

