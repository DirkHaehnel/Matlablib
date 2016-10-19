cd 'c:\Joerg\Doc\Nanophotonics&Optics\FrankSchleifenbaum\Pi Ethylenglycol\Lifetimes\'

clear all
close all

nglass = 1.52; % glass
nglycol = 1.43;
dagbottom = 0.036;
dagtop = 0.072;
dv = 0.025:0.001:0.25;

if 0 % Lifetime fitting
    load 'c:\Joerg\Doc\Nanophotonics&Optics\FrankSchleifenbaum\Pi Ethylenglycol\IRF.mat'

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

    ind = irf(:,1)>35 & irf(:,1)<100;
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

if 0 % Cavity transmission fitting
    cd 'c:\Joerg\Doc\Nanophotonics&Optics\FrankSchleifenbaum\Pi Ethylenglycol\Transmisssion\'

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

    para = [188 41];
    for j=1:size(y,3)
        para = Simplex('TransmissionFit',para,[],[],[],[],y(:,1,j),y(:,2,j));
        para = Simplex('TransmissionFit',para,[],[],[],[],y(:,1,j),y(:,2,j));
        res(j).d = para(1);
        res(j).dag = para(2);
        
        eval(['print -dpng -r300 ' fnames(j).name(1:end-4)])        
    end
    save TransmissionData res fnames
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
    load metals
    lamemv = 0.5:0.025:0.725;
    lamex = 0.47;

    nagex = interp1(SilverDmitri(:,1),SilverDmitri(:,2),1e3*lamex);
    nauex = interp1(wavelength,gold,1e3*lamex);
    rhofield = [0 1];
    NA = 1.25;
    fd = 164.5e3/100;
    over = 3e3;
    focpos = 0;
    mag = 100;
    av = 75;
    zpin = 0;
    atf = [];
    t = 0.01:0.01:20;
    n0 = [nglass nagex nauex]; d0 = [dagbottom 0.004];
    n = nglycol;
    n1 = [nauex nagex nglass]; d1 = [0.004 dagtop];
    for jd = 1:length(dv)
        d = dv(jd);
        resolution = [20 lamex/(d/100)];
        exc(jd) = GaussExc(rhofield, [0, d], NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, atf, resolution);
    end
    for jlamem=1:length(lamemv)
        lamem = lamemv(jlamem);
        nagem = interp1(SilverDmitri(:,1),SilverDmitri(:,2),1e3*lamem);
        nauem = interp1(wavelength,gold,1e3*lamem);

        for jd = 1:length(dv)
            d = dv(jd)/lamem*2*pi;

            n0 = [nglass nagem nauem];
            n1 = [nauem nagem nglass];

            mdf0(jd) = GaussExc2MDF(exc(jd), NA, n0, n, n1, focpos, lamem, mag, av, zpin, atf, 0, 3);
            mdf1(jd) = GaussExc2MDF(exc(jd), NA, n0, n, n1, focpos, lamem, mag, av, zpin, atf, 1, 3);
            [bla bla bla bla bla bla bla bla qv(:,jd),qp(:,jd)] = LifetimeL(exc(jd).z(1,:)/lamem*2*pi,n0,n,n1,d0/lamem*2*pi,d,d1/lamem*2*pi);
            jd
        end
        eval(['save QY_PerrylenimidGlycol' mint2str(lamem*1e3,3) ' exc mdf0 mdf1 qv qp'])
    end
end

if 0 % Cavity transmission
    load SilverDmitri
    load metals
    lamexv = 450:800;

    clear tp rp ts rs
    for jd = 1:length(dv)
        d = dv(jd);
        for jl = 1:length(lamexv)
            nagex = interp1(SilverDmitri(:,1),SilverDmitri(:,2),lamexv(jl));
            nauex = interp1(wavelength,gold,lamexv(jl));
            n0 = [nglass nagex nauex]; d0 = [dagbottom 0.004];
            n1 = [nauex nagex nglass]; d1 = [0.004 dagtop];
            [rp(jl,jd), rs(jl,jd), tp(jl,jd), ts(jl,jd)] = Fresnel(nglass,[n0 nglycol n1],1e3*[d0 d d1]/lamexv(jl)*2*pi);
        end
        jd
    end
    save TransmissionModel dv lamexv rp rs tp ts
    pcolor(dv*1e3,lamexv,abs(tp).^2); shading interp    
    xlabel('cavity thickness (nm)'); ylabel('wavelength (nm)');
end

if 0 % data fitting
    load PerrylenimidData
    load QY_PerrylenimidGlycol
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
    tauav(19) = [];
    

%     for j=1:length(res)
%         hell(j) = sum(res(j).A)/meastime(j);         
%     end
% 
%     p = Simplex('AffineFit',1,0,[],[],[],thickness', hell', 1e3*dv', brightness0',[],1);
%     [err, c, z] = AffineFit(p,thickness', hell', 1e3*dv, brightness0',[],1);
%     plot(thickness,z,thickness,hell,'o')
%     xlabel('cavity thickness (nm)'); ylabel('fluorescence intensity (a.u.)');

    p = Simplex('AffineFit',1,0,[],[],[],thickness', 1./tauav', 1e3*dv', 1./tau1',[],1);
    [err, c, z] = AffineFit(p,thickness', 1./tauav', 1e3*dv, 1./tau1',[],1);
    plot(thickness,1./z,thickness,tauav,'o')
    xlabel('cavity thickness (nm)'); ylabel('average lifetime (ns)');

%     p = Simplex('AffineFit',1,0,[],[],[],thickness', 1./tauav', 1e3*dv', 1./tau0',[],1);
%     [err, c, z] = AffineFit(p,thickness', 1./tauav', 1e3*dv, 1./tau0',[],1);
%     plot(thickness,1./z,thickness,tauav,'o')
%     xlabel('cavity thickness (nm)'); ylabel('average lifetime (ns)');
    
    quantumyield = c(2)/sum(c);
    lifetime_water = ntoluol/1.33/sum(c);
end

return

load QY_PerrylenimidGlycol625.mat
plot(1e3*dv,tau1,1e3*dv,tau0)
load QY_PerrylenimidGlycol575.mat
hold on
plot(1e3*dv,tau1,'--',1e3*dv,tau0,'--')
load QY_PerrylenimidGlycol535.mat
plot(1e3*dv,tau1,':',1e3*dv,tau0,':')
xlabel('cavity height (nm)'); ylabel('average lifetime (\tau_0)')
legend({'slow rot. 625 nm','fast rot. 625 nm','slow rot. 575 nm','fast rot. 575 nm','slow rot. 535 nm','fast rot. 535 nm'})
print -dpng -r300 LifetimeModel
hold off

return

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
