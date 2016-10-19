% clear all

if 1
    load metals
    lamex = 0.373;
    lamem = 0.450;

    ng = 1.52;
    n = 1.5;
    nauex = silver(wavelength==1e3*lamex);
    nauem = silver(wavelength==1e3*lamem);
    d0 = 0.03;
    dv = 0.05:0.001:0.25;
    d1 = 0.06;
    rhofield = [0 0.4];
    NA = 1.4;
    fd = 1.65e3;
    over = 3e3;
    focpos = 0.1;
    mag = 100;
    av = 50;
    zpin = 0;
    atf = [];
    kappa = 1;
    resolution = 50;
    delta = 1e-3;
	resolution = [20 round(lamex/(d/100))];

    for j=1:length(dv)
        d = dv(j);
        exc = GaussExc(rhofield, [0 d], NA, fd, [ng nauex], n, [nauex ng], d0, d, d1, lamex, over, focpos, atf, resolution);
        if mod(1e3*d,10)==0
            FocusImage2D(exc.rho*1e3,exc.z*1e3,cat(3,exc.fxc,exc.fxs));
            eval(['print -dpng -r300 ''FrankCavityIntensity' mint2str(d*1e3,3) ' nm'''])
        end
        int(j) = sum(exc.rho(:,1)'*(abs(exc.fxc(:,:,1)).^2 + abs(exc.fyc(:,:,1)).^2 + abs(exc.fzc(:,:,1)).^2));
        for k =1:2
            int(j) = int(j) + 0.5*sum(exc.rho(:,1)'*(abs(exc.fxc(:,:,k+1)).^2 + abs(exc.fyc(:,:,k+1)).^2 + abs(exc.fzc(:,:,k+1)).^2));
            int(j) = int(j) + 0.5*sum(exc.rho(:,1)'*(abs(exc.fxs(:,:,k)).^2 + abs(exc.fys(:,:,k)).^2 + abs(exc.fzs(:,:,k)).^2));
        end
    end

    %     for j=1:size(exc.z,2)
    %         [lvd,lvu,lpd,lpu,qvd,qvu,qpd,qpu,qv,qp] = LifetimeL(exc.z(1,j)/lamem*2*pi, [ng nauex], n, [nauex ng], d0/lamem*2*pi, d/lamem*2*pi, d1/lamem*2*pi);
    %         taup(j) = 4/3*n/(qpu+qpd);
    %         tauv(j) = 4/3*n/(qvu+qvd);
    %     end
    %     q=1-1e-2; % quantum yield of photochemsitry
    %     plot(exc.z(1,:)*1e3,q./(q+(1-q)./taup),exc.z(1,:)*1e3,q./(q+(1-q)./tauv));
    %     xlabel('cavity height (nm)'); ylabel('quantum yield');
    
    plot(dv*1e3,int); xlabel('cavity height (nm)'); ylabel('total light intensity (a.u.)');
    
end
