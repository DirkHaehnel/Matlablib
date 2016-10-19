if 0
    names = {'atto655_10uW.pt3',...
        'atto655_50uW.pt3',...
        'atto655_100uW.pt3',...
        'atto655_200uW.pt3',...
        'atto655_300uW.pt3',...
        'atto655_400uW.pt3'};
    
    pow = [10.0912 49.8848 100.2456 201.4432 302.1172 404.7428];

    dname = 'c:\Joerg\Doc\Fcs\Photophysics\Atto655\';

    for j=1:length(names)
        [y, tcspc, chan, markers, num, overcount, head] = pt3v2read([dname names{j}], [1 inf]);
        x(:,j) = mhist(tcspc*head.Resolution,25:25:2e3);
        t(j) = head.StopAfter/1e3;
    end
    [y, tcspc, chan, markers, num, overcount, head] = pt3v2read([dname 'backscattering.pt3'], [1 inf]);
    [irf, tau] = mhist(tcspc*head.Resolution,25:25:2e3);


    plot(tau(2:end-1),x(2:end-1,:)./(ones(size(x,1)-2,1)*max(x)),tau(2:end-1),irf(2:end-1,:)./max(irf),'o')
    ylabel('fluorescence intensity (a.u.)')
    xlabel('time (ns)')
    
    %y = x./(irf*ones(1,size(x,2)));
    %plot(tau(2:end-1),y(2:end-1,:)./(ones(size(y,1)-2,1)*max(y)))

    save AttoSatData x irf tau names pow t
    
    plot(pow,max(x).*[1 1 1 2 2 2],'-o')
    xlabel('excitation intensity (µW)')
    ylabel('max. fluorescence intensity (a.u.)')
    grid
    
    load c:\Joerg\Doc\Fcs\Photophysics\Atto655\Point_1.mat
    tmp = sum(res.tcspc);
    kappa = (sqrt(tmp(1)*tmp(4))-sqrt(tmp(2)*tmp(3)))/(sqrt(tmp(1)*tmp(4))+2*sqrt(tmp(2)*tmp(3)));
    for k=1:4
        [ind ind] = max(res.tcspc(:,k));
        close; 
        tau(k) = Simplex('ExpFun',1,0,[],[],[],res.tau(ind+50:ind+500,k),res.tcspc(ind+50:ind+500,k),1);
    end
end

if 0
    %x=10.^(-4:0.1:-0.5);
    %plot(x,PulsedExcitation(x,1.84,0.05,25),x,1/25*(1-exp(-x*25)),x,1/25*ones(size(x)),':');
    clear all
    
    NA = 1.2;
    fd = 3e3;
    n0 = 1.333;
    n = 1.333;
    n1 = 1.333;
    d0 = [];
    d = [];
    d1 = [];
    lamex = 0.64;
    overv = (1.1:0.05:1.4)*1e3;
    focpos = 20;
    av = 150/2; % confocal pinhole radius
    lamem = 0.670;
    mag = 60;
    zpin = 0;
    atf = [];
    resolution = 2*[20 20];
    kappa = 0.06/(2/5); % anisotropy factor
    lt = [];
    lifetime = 1.84; % lifetime
    pulse = [0.05 25]/lifetime;  % pulse width and rep. period
    satv = 10.^(-4:0.1:1);
    triplet = [];

    rhofield = [0 20];
    zfield = [0 40];

    for jo=1:length(overv)
        over = overv(jo);
        exc = GaussExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, atf, resolution);
        wave = abs(exc.fxc(:,:,1)).^2 + abs(exc.fyc(:,:,1)).^2 + abs(exc.fzc(:,:,1)).^2 + sum(abs(exc.fxc(:,:,2:end)).^2 + abs(exc.fyc(:,:,2:end)).^2 + abs(exc.fzc(:,:,2:end)).^2,3)/2;
        intx = 0*satv; inty = 0*satv;
        for k=1:length(satv)
            sat = satv(k);
            disp(sat)
            mdf = GaussExc2MDF(exc, NA, n0, n, n1, focpos, lamem, mag, av, zpin, atf, kappa, lt, pulse, sat, triplet);
            intx(k) = sum(exc.rho(:,1)'*mdf.volx(:,:,1));
            inty(k) = sum(exc.rho(:,1)'*mdf.voly(:,:,1));

            %FocusImage3D(exc.rho,exc.z-focpos,mdf.volx+mdf.voly,[],[],3*pi/2)
            %view([25 30])
            %eval(['print -dpng -r300 AttoSatMDF' mint2str(k,2)])
        end
        s0 = 2*pi*max(exc.rho(:,1)'*(abs(exc.fxc(:,:,1)).^2 + abs(exc.fyc(:,:,1)).^2 + abs(exc.fzc(:,:,1)).^2))*diff(exc.rho(1:2,1));
        eval(['save AttoSatTheory' mint2str(jo,1)])
    end
end
 
if 0 % data fitting
    %load c:\Joerg\MATLAB\AttoSatData.mat pow x
    tmp = load('c:\Joerg\Doc\Fcs\DSOM-FCS\AttoSaturation.txt');
    pow = tmp(:,1)'; x = tmp(:,2)';
    load c:\Joerg\MATLAB\AttoSatTheory.mat satv intx inty exc s0
    load c:\Joerg\Doc\Fcs\Photophysics\Atto655\attoabsorption.mat
    NatConst
    
    sigma = 3*1e11*log(10)*attoabs(attoabs(:,1)==640,2)*extinction/AvogadroConstant; % sigma in mum^2
    intsat = PlanckConstant*SpeedOfLight/640e-9/sigma/(lifetime*1e-9)*1e6; % saturation intensity in muW/mum^2
    peak_enhancement = max(abs(exc.fxc(1,:,1)).^2 + abs(exc.fyc(1,:,1)).^2 + abs(exc.fzc(1,:,1)).^2)/s0;
    
    close;
    int = satv.*(intx+inty);

    k = length(pow);
    para = Simplex('AffineFit',1,0,[],[],[],[-pow(k:-1:1)';0;pow(1:k)'],[-flipud(x(1:k)');0;x(1:k)'],...
        intsat/peak_enhancement*[-satv(end:-1:1)';satv'],[-flipud(int');int'],[],2);
    
    [err, c, z, zz] = AffineFit(para,[-pow(k:-1:1)';0;pow(1:k)'],[-flipud(x(1:k)');0;x(1:k)'],...
        intsat/peak_enhancement*[-satv(end:-1:1)';satv'],[-flipud(int');int'],[],2);

    %plot(para*intsat/peak_enhancement*satv,c(2)*int,pow,max(x).*[1 1 1 2 2 2],'o')
    close
    plot(pow,x,'ob')
    xlabel('excitation intensity (µW)')
    ylabel('max. fluorescence intensity (a.u.)')
    ax1 = gca;
    ax2 = axes('Position',get(ax1,'Position'),...
        'XAxisLocation','top',...
        'Color','none',...
        'XColor','k','YColor','k');
    line(para*intsat/peak_enhancement*satv,c(2)*int)
    set(ax2,'xlim',get(ax1,'xlim'))
    tmp = get(ax1,'xtick')/(para*intsat/peak_enhancement);
    clear ss; ss{1}='0';
    for j=2:length(tmp)
        ss{j} = mnum2str(tmp(j),1,2);
    end
    set(ax2,'xticklabel',ss)
    grid
    
end



