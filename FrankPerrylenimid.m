cd 'c:\Joerg\Doc\Nanophotonics&Optics\FrankSchleifenbaum\Liquid Cavity\Lifetimes'

clear all
close all

nglass = 1.52; % glass
ntoluol = 1.496;
dagbottom = 0.036;
dagtop = 0.072;
dv = 0.025:0.001:0.25;

if 0 % Lifetime fitting
    load c:\Joerg\Doc\Nanophotonics&Optics\FrankSchleifenbaum\EYFP\irf.mat % irf from 6.11.07

    fnames = dir('*.txt');
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
    tau = [0.5 2 5];
    irf = circshift(irf(ind,2),0);
    for j=1:size(y,3)
        close all

        [c, offset, A, tau] = fluofit(irf, y(ind,2,j), p, dt);
        [c, offset, A, tau] = fluofit(irf, y(ind,2,j), p, dt, [0.001 tau']);
        [c, offset, A, tau] = fluofit(irf, y(ind,2,j), p, dt, tau');
        res(j).c = c;
        res(j).offset = offset;
        res(j).A = A;
        res(j).tau = tau;

        eval(['print -dpng -r300 ' fnames(j).name(1:end-4)])
    end
    save PerrylenimidData res fnames
end

if 0 % average lifetime estimating
    load c:\Joerg\Doc\Nanophotonics&Optics\FrankSchleifenbaum\EYFP\irf.mat % irf from 6.11.07
    
    fnames = dir('*.txt');
    for j=1:length(fnames)
        fin = fopen(fnames(j).name,'r');
        for k=1:3
            fgetl(fin);
        end
        tmp = fscanf(fin,'%f %f\n',inf);
        fclose(fin);
        y(:,:,j) = [tmp(1:2:end) tmp(2:2:end)];
    end

    ind = irf(:,1)>35 & irf(:,1)<75;
    tp = irf(ind,1);
    % semilogy(irf(ind,1),irf(ind,2),irf(ind,1),ones(sum(ind),1)*mean(irf(irf(:,1)>35 & irf(:,1)<37,2)))
    % for j=1:51 semilogy(y(ind,1,j),y(ind,2,j),y(ind,1,j),ones(sum(ind),1)*mean(y(irf(:,1)>35 & irf(:,1)<37,2,j))); pause; end
    zz = irf(ind,2)-mean(irf(irf(:,1)>35 & irf(:,1)<37,2));
    h1 = sum(zz.*tp);
    h0 = sum(zz);
    clear tauav
    for j=1:size(y,3)
        yy = y(ind,2,j)-mean(y(irf(:,1)>35 & irf(:,1)<37,2,j));
        f1 = sum(yy.*tp);
        f0 = sum(yy);
        tauav(j) = f1/f0-h1/h0;
    end
    save PerrylenimidData tauav -append
end

if 1 % Cavity lifetime modeling
    load SilverDmitri
    lamem = 0.535;
    lamex = 0.47;

    nagex = interp1(SilverDmitri(:,1),SilverDmitri(:,2),1e3*lamex);
    nagem = interp1(SilverDmitri(:,1),SilverDmitri(:,2),1e3*lamem);

    rhofield = [0 1];
    NA = 1.25;
    fd = 164.5e3/100;
    over = 3e3;
    focpos = 0;
    mag = 100;
    av = 75;
    zpin = 0;
    atf = [];
    kappa = 1;
    t = 0.01:0.01:20;
    for jd = 1:length(dv)
        d = dv(jd);

        n0 = [nglass nagex]; d0 = dagbottom;
        n = ntoluol; 
        n1 = [nagex nglass]; d1 = dagtop;
        
        resolution = [20 lamex/(d/100)];
        exc = GaussExc(rhofield, [0, d], NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, atf, resolution);
        
        n0 = [nglass nagem];
        n1 = [nagem nglass];
        
        mdf = GaussExc2MDF(exc, NA, n0, n, n1, focpos, lamem, mag, av, zpin, atf, kappa, 1);
        tau = GaussExc2MDF(exc, NA, n0, n, n1, focpos, lamem, mag, av, zpin, atf, kappa, 2);
        tau1(jd) = 4/3*n*squeeze(sum(sum(exc.rho.*(tau.volx(:,:,1)+tau.voly(:,:,1)))))/squeeze(sum(sum(exc.rho.*(mdf.volx(:,:,1)+mdf.voly(:,:,1)))));
        brightness1(jd) = squeeze(sum(sum(exc.rho.*(mdf.volx(:,:,1)+mdf.voly(:,:,1)))))*diff(exc.z(1,1:2));

        mdf = GaussExc2MDF(exc, NA, n0, n, n1, focpos, lamem, mag, av, zpin, atf, 0, 1);
        %[lvd(jd,:),lvu(jd,:),lpd(jd,:),lpu(jd,:),qvd(jd,:),qvu(jd,:),qpd(jd,:),qpu(jd,:)] = LifetimeL(2*pi*exc.z(1,:)/lamem, n0, n, n1, 2*pi*d0/lamem, 2*pi*d/lamem, 2*pi*d1/lamem);
        %lifetime(:,jd) = 4*n./(2*(qpu(jd,:)+qpd(jd,:))+qvu(jd,:)+qvd(jd,:))';
        qv(jd,:) = mdf.qv; qp(jd,:) = mdf.qp;
        lifetime(:,jd) = 4*n./(2*mdf.qp+mdf.qv)';
        weight(:,jd) = squeeze(sum(exc.rho.*(mdf.volx(:,:,1)+mdf.voly(:,:,1))))';
        tau0(jd) = weight(:,jd)'*lifetime(:,jd)/sum(weight(:,jd));
        brightness0(jd) = squeeze(sum(sum(exc.rho.*(mdf.volx(:,:,1)+mdf.voly(:,:,1)))))*diff(exc.z(1,1:2));
        jd
    end
    save QY_Perrylenimid lifetime weight tau0 brightness0 qv qp
    %save QY_Perrylenimid lifetime weight tau0 tau1 brightness0 brightness1 lvd lvu lpd lpu qvd qvu qpd qpu
    % pcolor(ones(size(lifetime,1),1)*dv,(0.1:size(lifetime,1))'*dv,lifetime); shading flat
    % pcolor(ones(size(lifetime,1),1)*dv,(0.1:size(lifetime,1))'*dv,4/3*ntoluol./(qvu+qvd)'); shading flat
    % pcolor(ones(size(lifetime,1),1)*dv,(0.1:size(lifetime,1))'*dv,4/3*ntoluol./(qpu+qpd)'); shading flat
    % pcolor(ones(size(lifetime,1),1)*dv,(0.1:size(lifetime,1))'*dv,weight); shading flat
end

if 0 % Cavity transmission
    load SilverDmitri
    lamexv = 450:850;

    clear tp rp ts rs
    for jd = 1:length(dv)
        d = dv(jd);
        for jl = 1:length(lamexv)
            nagex = interp1(SilverDmitri(:,1),SilverDmitri(:,2),lamexv(jl));
            n0 = [nglass nagex]; d0 = dagbottom;
            n1 = [nagex nglass]; d1 = dagtop;
            [rp(jl,jd), rs(jl,jd), tp(jl,jd), ts(jl,jd)] = Fresnel(nglass,[n0 ntoluol n1],1e3*[d0 d d1]/lamexv(jl)*2*pi);
        end
        jd
    end
    save TransmissionModel dv lamexv rp rs tp ts
    pcolor(dv*1e3,lamexv,abs(tp).^2); shading interp    
    xlabel('cavity thickness (nm)'); ylabel('wavelength (nm)');
end

if 0 % data fitting
    load PerrylenimidData
    load QY_Perrylenimid
    load TransmissionModel
    load TransmissionData

    p = 100;
    for jd=1:size(tp,2)
        maxlam(jd) = lamexv(abs(tp(:,jd))==max(abs(tp(:,jd))));
    end
    ind = dv>0.1 & dv<0.2;
    ct = polyfit(maxlam(ind),1e3*(dv(ind)),1);
    thickness = polyval(ct,transmission);
    
	for j=1:length(res) 
        ind = res(j).tau<10.0; 
        w = res(j).A(ind)/sum(res(j).A(ind));
        %tauav(j)=1/(w'*(1./res(j).tau(ind))); 
        tauav(j) = w'*res(j).tau(ind); 
    end

    for j=1:length(res)
        hell(j) = sum(res(j).A)/meastime(j);         
    end

    p = Simplex('AffineFit',1,0,[],[],[],thickness', hell', 1e3*dv', brightness0',[],1);
    [err, c, z] = AffineFit(p,thickness', hell', 1e3*dv, brightness0',[],1);
    plot(thickness,z,thickness,hell,'o')
    xlabel('cavity thickness (nm)'); ylabel('fluorescence intensity (a.u.)');

    p = Simplex('AffineFit',1,0,[],[],[],thickness', 1./tauav', 1e3*dv', 1./tau1',[],1);
    [err, c, z] = AffineFit(p,thickness', 1./tauav', 1e3*dv, 1./tau1',[],1);
    plot(thickness,1./z,thickness,tauav,'o')
    xlabel('cavity thickness (nm)'); ylabel('average lifetime (ns)');

    p = Simplex('AffineFit',1,0,[],[],[],thickness', 1./tauav', 1e3*dv', 1./tau0',[],1);
    [err, c, z] = AffineFit(p,thickness', 1./tauav', 1e3*dv, 1./tau0',[],1);
    plot(thickness,1./z,thickness,tauav,'o')
    xlabel('cavity thickness (nm)'); ylabel('average lifetime (ns)');
    
    quantumyield = c(2)/sum(c);
    lifetime_water = ntoluol/1.33/sum(c);
end

