% program for GFP saturation modeling

if 0 % membrane
    rhofield = [0 1];
    zfield = 0;
    NA = 1.2;
    fd = 3e3;
    n0 = 1.33;
    n = n0;
    n1 = n0;
    d0 = [];
    d = 0;
    d1 = [];
    lamex = 0.488;
    over = 3e3;
    focpos = 0;
    atf = [];
    resolution = 50;

    lamem = 0.53;
    mag = 60;
    av = 50;
    zpin = 0;

    exc = GaussExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, atf, resolution);

    phi=0:pi/50:2*pi;
    ff = abs(exc.fxc(:,:,1)*cos(0*phi)+exc.fxc(:,:,3)*cos(2*phi)).^2 + abs(exc.fys(:,:,2)*sin(2*phi)).^2 + abs(exc.fzc(:,:,2)*cos(phi)).^2;
    ff = ff/max(ff(:));

    mdf = GaussExc2MDF(exc, NA, n0, n, n1, focpos, lamem, mag, av, zpin);

    gg = (mdf.volx(:,:,1) + mdf.voly(:,:,1))*cos(0*phi);
    for j=1:(size(mdf.volx,3)-1)/2
        gg = gg + (mdf.volx(:,:,j+1) + mdf.voly(:,:,j+1))*cos(j*phi);
        gg = gg + (mdf.volx(:,:,(end+1)/2+j) + mdf.voly(:,:,(end+1)/2+j))*cos(j*phi);
    end
    gg = gg./ff; gg = gg/max(gg(:));

    k10 = 1;
    Tpulse = 0.05/3
    Trep = 25/3;
    satv = 0:0.01:1;
    for j=1:length(satv)
        y = PulsedExcitation(ff*satv(j),k10,Tpulse,Trep)*Trep;
        for k=1:5
            int(j,k) = sum(exc.rho'*(1-(1-y).^k));
        end
    end

    plot(satv,int./(ones(length(satv),1)*diff(int(1:2,:))))
end


if 1 % solution
    rhofield = [0 6];
    zfield = [0 10];
    NA = 1.2;
    fd = 3e3;
    n0 = 1.33;
    n = n0;
    n1 = n0;
    d0 = [];
    d = 0;
    d1 = [];
    lamex = 0.488;
    over = 3e3;
    focpos = 0;
    atf = [];
    resolution = [50 10];

    lamem = 0.53;
    mag = 60;
    av = 50;
    zpin = 0;

    exc = GaussExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, atf, resolution);

    ff = abs(exc.fxc(:,:,1)).^2;
    ff = ff/max(ff(:));

    mdf = GaussExc2MDF(exc, NA, n0, n, n1, focpos, lamem, mag, av, zpin);

    gg = (mdf.volx(:,:,1) + mdf.voly(:,:,1));
    gg = gg./ff; gg = gg/max(gg(:));

    k10 = 1;
    Tpulse = 0.05/3
    Trep = 25/3;
    satv = 0:1:100;
    for j=1:length(satv)
        y = PulsedExcitation(ff*satv(j),k10,Tpulse,Trep)*Trep;
        for k=1:5
            int(j,k) = sum(exc.rho(:,1)'*(1-(1-y).^k));
        end
    end

    plot(satv,int./(ones(length(satv),1)*diff(int(1:2,:))))
end

