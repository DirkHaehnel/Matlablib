clear
close all

if 0 % PSF figure
    over = 1.25e3;
    cov = 0;
    abev = [];
    rhofield = [0 3];
    zfield = [-5 5];
    NA = 1.2;
    n0 = 1.333;
    n = 1.333;
    n1 = 1.333;
    d0 = [];
    d = 0;
    d1 = [];
    lamex = 0.64;
    focpos1 = 0;
    focpos2 = 0.3;
    lamem = 0.67;
    mag = 60;
    av = 75;
    zpin = 0e3;
    atf = [];
    lt = [];
    kappa = 1;
    fd = 3e3;
    resolution = 2*[20 20];

    exc = GaussExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos1, atf, resolution);
    pulse = [0.05/2 25/2];
    sat = 0.01;
    mdf1 = GaussExc2MDF(exc, NA, n0, n, n1, focpos1, lamem, mag, av, zpin, atf, kappa, lt, pulse, sat);
    sat = 0.1;
    mdf2 = GaussExc2MDF(exc, NA, n0, n, n1, focpos2, lamem, mag, av, zpin, atf, kappa, lt, pulse, sat);
    rho = exc1.rho;
    z = exc1.z;

    figure
    s = 1;
    subplot(1,3,1)
    tmp1 = (mdf1.volx+mdf1.voly);
    FocusImage2D(rho,z,tmp1);
    subplot(1,3,2)
    tmp2 = (mdf2.volx+mdf2.voly);
    FocusImage2D(rho,z,tmp2);
    subplot(1,3,3)
    tmp12 = tmp1-s*tmp2;
    FocusImage2D(rho,z,tmp12);
end

if 0 % modeling
    over = 1.25e3;
    cov = 0;
    abev = [];
    rhofield = [0 3];
    zfield = [-5 5];
    NA = 1.2;
    n0 = 1.333;
    n = 1.333;
    n1 = 1.333;
    d0 = [];
    d = 0;
    d1 = [];
    lamex = 0.64;
    focpos1 = 0;
    focpos2 = 0;
    lamem = 0.67;
    mag = 60;
    av = 75;
    zpin = 0e3;
    atf = [];
    lt = [];
    kappa = 1;
    fd = 3e3;
    resolution = 2*[20 20];
    maxm = 10;

    exc = GaussExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos1, atf, resolution, [], maxm);
    pulse = [0.05/2 25/2];
    sat = 0.;
    mdf0 = GaussExc2MDF(exc, NA, n0, n, n1, focpos1, lamem, mag, av, zpin, atf, kappa, lt, pulse, sat);
    [detvol0, intens0] = DetectionVolume(exc.rho,exc.z,mdf0.volx(:,:,1),mdf0.voly(:,:,1));
    [modres0, autotime] = FCS(exc.rho,exc.z,mdf0.volx,mdf0.voly);
    tau0 = exp(interp1(modres0(:,1),log(autotime),0.5,'cubic'));
    satv = 0.01*exp((0:20)/20*log(0.3/0.01));
    %exc = GaussExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos2, atf, resolution);
    exc = GaussExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, [focpos1 0.005], atf, resolution, [], maxm);
    clear mdf
    for k=1:length(satv)
        sat = satv(k);
        mdf(k) = GaussExc2MDF(exc, NA, n0, n, n1, focpos2, lamem, mag, av, zpin, atf, kappa, lt, pulse, sat);
        [modres(:,:,k), autotime] = FCS(exc.rho,exc.z,mdf0.volx,mdf0.voly,mdf(k).volx,mdf(k).voly);
        sv = 0:0.005:2;
        for j=1:length(sv)
            s = sv(j);
            modres12(:,j,k) = modres(:,1,k) - 2*s*modres(:,3,k) + s^2*modres(:,2,k);
            modres12(:,j,k) = modres12(:,j,k)/max(modres12(:,j,k));
            [detvol(j,k), intens(j,:,k)] = DetectionVolume(exc.rho,exc.z,mdf0.volx(:,:,1)-s*mdf(k).volx(:,:,1),mdf0.voly(:,:,1)-s*mdf(k).voly(:,:,1));
        end

        for j=1:length(sv)
            tau2(j,k) = exp(interp1(modres12(:,j,k),log(autotime),0.5,'cubic'));
        end
    end

    save OptSatMicroscopeX005
end

if 0 % figures
    % diffusion time plot
    plot(sv,tau2/tau0)
    ax=axis;
    axis([sv([1 end]) ax(3:4)]);
    xlabel('subtracting factor'); ylabel('rel. diffusion time \tau_{\itD\rm,1/2}')
    % print -dpng -r300 OptSatMicroscopeTau

    % intensity plot
    plot(sv,mean(intens,2)/mean(intens0))
    ax=axis;
    axis([sv([1 end]) ax(3:4)]);
    grid
    xlabel('subtracting factor'); ylabel('rel. fluorescence intensity')
    % print -dpng -r300 OptSatMicroscopeIntensity

    % detection volume
    plot(sv,detvol/detvol0)
    ax=axis;
    axis([sv([1 end]) ax(3:4)]);
    grid
    xlabel('subtracting factor'); ylabel('rel. detection volume')
    % print -dpng -r300 OptSatMicroscopeDeteVol

    % ACFs
    surf(log10(autotime)*ones(1,length(sv)),ones(size(modres12,1),1)*sv,modres12,'facealpha',0.5)
    shading interp
    hold on; mesh(log10(autotime)*ones(1,length(sv(1:10:end))),ones(size(modres12,1),1)*sv(1:10:end),modres12(:,1:10:end),'facealpha',0.5); hold off
    %set(gca,'xscale','log')
    axis tight
    xlabel('time (a.u.)'); ylabel('subtracting factor'); zlabel('autocorrelation')
    cameratoolbar
    camlight
    % print -dpng -r300 OptSatMicroscopeACF
end

if 1 % data fitting
    load('c:\Joerg\Doc\Fcs\DSOMFCS\DSOMFCSdata.mat')
    svdata = sv; taudata = tau;
    load OptSatMicroscope 

    shiftv = 0:1e-4:5e-3;
    err = 0*shiftv;
    for j=1:length(shiftv)
        pdet = Simplex('AffineFit',0.9,0,[],[],[],data(ind,1)-xd-shiftv(j),data(ind,2),sv-xm,detvol/detvol0,[],2);
        err(j) = AffineFit(pdet,data(ind,1)-xd-shiftv(j),data(ind,2),sv-xm,detvol/detvol0,[],2);
    end
    shift = shiftv(err==min(err));
    pdet = Simplex('AffineFit',0.9,0,[],[],[],data(ind,1)-xd-shift,data(ind,2),sv-xm,detvol/detvol0,[],2);
    [err, c] = AffineFit(pdet,data(ind,1)-xd-shift,data(ind,2),sv-xm,detvol/detvol0,[],2);
    plot(sv,detvol/detvol0,xm+(data(ind,1)-xd-shift)/pdet,data(ind,2)/c(2),'o');
    axis([xm+(data(ind([1 end]),1)'-xd-shift)/pdet 0 1])
    xlabel('subtracting factor');
    ylabel('rel. detection volume');

    [err, c] = AffineFit(pdet,data(ind,1)-xd-shift,data(ind,4),sv-xm,tau2/tau0,[],2);
    plot(sv,c(1)+c(2)*tau2/tau0,xm+(data(ind,1)-xd-shift)/pdet,data(ind,4),'o');
    axis([xm+(data(ind([1 end]),1)'-xd-shift)/pdet 0 0.4])
    xlabel('subtracting factor');
    ylabel('diffusion time (ms)');

    [err, c] = AffineFit(pdet,data(ind,1)-xd-shift,data(ind,6),sv-xm,mean(intens,2)/mean(intens0),[],2);
    plot(sv,c(1)+c(2)*mean(intens,2)/mean(intens0),xm+(data(ind,1)-xd-shift)/pdet,data(ind,6),'o');
    axis([xm+(data(ind([1 end]),1)'-xd-shift)/pdet -3e3 4e3])
    xlabel('subtracting factor');
    ylabel('diffusion time (ms)');
end