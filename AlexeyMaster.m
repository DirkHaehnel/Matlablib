if 0
    [tau, lam, data] = AlexeyLifetimeRead('c:\Joerg\Doc\Nanophotonics&Optics\AlexeyChizhik\RotoDiffCavity\Laser signal.txt','c:\Joerg\Doc\Nanophotonics&Optics\AlexeyChizhik\RotoDiffCavity\Oregon green\');
    save 'c:\Joerg\Doc\Nanophotonics&Optics\AlexeyChizhik\RotoDiffCavity\LifetimeOregongreen.mat' tau lam data
    
    tau0 = AlexeyLifetimeRead('c:\Joerg\Doc\Nanophotonics&Optics\AlexeyChizhik\RotoDiffCavity\Laser signal.txt','c:\Joerg\Doc\Nanophotonics&Optics\AlexeyChizhik\RotoDiffCavity\free space lifetime oregongreen.txt');
    save 'c:\Joerg\Doc\Nanophotonics&Optics\AlexeyChizhik\RotoDiffCavity\LifetimeOregongreen.mat' tau0 -append
end

if 0
    dv = 100:5:250;
    for j=40:40
        model  = AlexeyModel('silver', 1.33, 488, j, 85, dv); eval(['save ModelSilver488ex' mint2str(j,2) 'nm model']);
    end
    for j=30:5:40
        model  = AlexeyModel('silver', 1.47, 488, j, 85, dv); eval(['save ModelSilver488exGlycerol' mint2str(j,2) 'nm model']);
    end
    for j=30:5:40
        model  = AlexeyModel('gold', 1.33, 488, j, 85, dv); eval(['save ModelGold488ex' mint2str(j,2) 'nm model']);
    end
    for j=30:5:40
        model  = AlexeyModel('gold', 1.47, 488, j, 85, dv); eval(['save ModelGold488exGlycerol' mint2str(j,2) 'nm model']);
    end
end


if 0
    load('c:\Joerg\Doc\Nanophotonics&Optics\AlexeyChizhik\RotoDiffCavity\LifetimeOregongreen.mat');
    ltdata = [lam(:) tau(:)];
    ltdata(end,:) = []; % last data point seems to be wrong
    spectrum = load('c:\Joerg\Doc\Nanophotonics&Optics\AlexeyChizhik\RotoDiffCavity\free space spectrum oregongreen.txt');
    spectrum = spectrum(spectrum(:,1)>500,:); % bandpass filter
    load('c:\Joerg\Doc\Nanophotonics&Optics\AlexeyChizhik\RotoDiffCavity\ModelSilver488ex30nm.mat','model');
    
    [para, facv, ind] = AlexeyCavityFit(ltdata, spectrum, model, tau0);
end

if 0
    load('c:\Joerg\Doc\Nanophotonics&Optics\AlexeyChizhik\RotoDiffCavity\LifetimeFluorescein.mat');
    ltdata = [lam(:) tau(:)];
    ltdata(end,:) = []; % last data point seems to be wrong
    spectrum = load('c:\Joerg\Doc\Nanophotonics&Optics\AlexeyChizhik\RotoDiffCavity\free space spectrum fluorescein.txt');
    spectrum = spectrum(spectrum(:,1)>500,:); % bandpass filter
    load('c:\Joerg\Doc\Nanophotonics&Optics\AlexeyChizhik\RotoDiffCavity\ModelSilver488ex30nm.mat','model');
    
    [para, facv, ind] = AlexeyCavityFit(ltdata, spectrum, model, tau0);
end

if 0
    ltdata = load('c:\Joerg\Doc\Nanophotonics&Optics\AlexeyChizhik\RotoDiffCavity\Water R6G.txt');
    
    %ltdata = ltdata(7:end-1,:);
    
    spectrum = load('c:\Joerg\Doc\Nanophotonics&Optics\AlexeyChizhik\RotoDiffCavity\spectrum R6G in water.txt');
    tau0 = 4;
    load('c:\Joerg\Doc\Nanophotonics&Optics\AlexeyChizhik\RotoDiffCavity\ModelSilver488ex30nm.mat','model');
    
    [para, facv, ind] = AlexeyCavityFit(ltdata, spectrum, model, tau0);
end

if 0 % glycerol
    nsolvent = 1.47;
    dagtop = 85;
    lamex = 488;
    
    % fitting the transmission curves
    dname = 'c:\Joerg\Doc\Nanophotonics&Optics\AlexeyChizhik\RotoDiffCavity\Converted transmissions\R6G in glycerol 20111121\';
    fnames = dir([dname '*.asc']);
    load('c:\Joerg\Doc\Nanophotonics&Optics\AlexeyChizhik\RotoDiffCavity\ModelSilver488exGlycerol30nm.mat','model');
    irf = load('c:\Joerg\Doc\Nanophotonics&Optics\AlexeyChizhik\RotoDiffCavity\TransmissionGauge.txt');
    load('c:\Joerg\Doc\Nanophotonics&Optics\AlexeyChizhik\RotoDiffCavity\ag_ellipsometry.mat');
    ind = 5:25;
    p = zeros(2,numel(fnames));
    for j=1:length(fnames)
        x = load([dname fnames(j).name]);
        x(:,2) = (x(:,2)-min(x(:,2)))./interp1(irf(:,1),irf(:,2),x(:,1));
        mm = x(x(:,2)==max(x(:,2)),1);
        close
        %p(:,j) = Simplex('CavityTransmissionFit', [polyval(polyfit(model.tmax(ind),model.dv(ind),1),mean(mm)) 30], [], [], [], [], x(:,1), x(:,2), nsolvent, 'silver', dagtop);
        p(:,j) = Simplex('CavityTransmissionFit', [polyval(polyfit(model.tmax(ind),model.dv(ind),1),mean(mm)) 30], [], [], [], [], x(:,1), x(:,2), nsolvent, [ag_lam ag_re+1i*ag_im], dagtop);
        num(j) = str2double(fnames(j).name(1:end-4));
        xlabel('wavelength (nm)'); ylabel('transmission (a.u.)');
        ax = axis;
        s = {['\itd\rm = ' mnum2str(p(1,j),3,1) ' nm'],['\itd\rm_{Ag} = ' mnum2str(p(2,j),3,1) ' nm']};
        text(ax(1)+0.7*diff(ax(1:2)),ax(3)+0.8*diff(ax(3:4)),s);
        eval(['print -dpng -r300 ' '''' dname fnames(j).name(1:end-4) '''' 'New'])
    end
    x = load([dname 'transmissions and decays.txt']);
    dv = zeros(size(x,1),2); dag = dv; tmax = zeros(size(x,1),1);
    for j=1:size(x,1)
        dv(j,:) = [p(1,num==x(j,1)) p(1,num==x(j,2))];
        dag(j,:) = [p(2,num==x(j,1)) p(2,num==x(j,2))];
        tmax(j) = x(j,3);
    end
	eval(['save ' '''' dname 'TransmissionFitNew' '''' ' tmax dv dag'])

    % modeling
    dv = mean(dv,2);
    dag = mean(dag,2);
    %model  = AlexeyModel('silver', nsolvent, lamex, dag, dagtop, dv);
    model  = AlexeyModel([ag_lam ag_re+1i*ag_im], nsolvent, lamex, dag, dagtop, dv);
    model.tmax = tmax;
    eval(['save ' '''' dname 'modelNew' '''' ' model'])
    
    ltdata = load('c:\Joerg\Doc\Nanophotonics&Optics\AlexeyChizhik\RotoDiffCavity\Glycerol R6G.txt');
    spectrum = load('c:\Joerg\Doc\Nanophotonics&Optics\AlexeyChizhik\RotoDiffCavity\spectrum R6G in glycerol.txt');
    tau0 = 3.9;
    load('c:\Joerg\Doc\Nanophotonics&Optics\AlexeyChizhik\RotoDiffCavity\Converted transmissions\R6G in glycerol 20111121\modelNew.mat','model');
    
    para = AlexeyCavityFitFull(ltdata, spectrum, model, []);
    print -dpng -r300 R6GGlycerolFitFullNew
    ind = ltdata(:,1)>=520 & ltdata(:,1)<=640;
    [parag taufitg lamg taug] = AlexeyCavityFitFull(ltdata(ind,:), spectrum, model, []);
    print -dpng -r300 R6GGlycerolFitReducedNew
end

if 1 % water
    nsolvent = 1.33;
    dagtop = 85;
    lamex = 488;

    % fitting the transmission curves
    dname = 'c:\Joerg\Doc\Nanophotonics&Optics\AlexeyChizhik\RotoDiffCavity\Converted transmissions\R6G in water 20111121\';
    fnames = dir([dname '*.asc']);
    load('c:\Joerg\Doc\Nanophotonics&Optics\AlexeyChizhik\RotoDiffCavity\ModelSilver488ex30nm.mat','model');
    irf = load('c:\Joerg\Doc\Nanophotonics&Optics\AlexeyChizhik\RotoDiffCavity\TransmissionGauge.txt');
    load('c:\Joerg\Doc\Nanophotonics&Optics\AlexeyChizhik\RotoDiffCavity\ag_ellipsometry.mat');
    ind = 5:25;
    p = zeros(2,numel(fnames));
    num = zeros(numel(fnames),1);
    for j=1:length(fnames)
        x = load([dname fnames(j).name]);
        x(:,2) = (x(:,2)-min(x(:,2)))./interp1(irf(:,1),irf(:,2),x(:,1));
        mm = x(x(:,2)==max(x(:,2)),1);
        close
        %p(:,j) = Simplex('CavityTransmissionFit', [polyval(polyfit(model.tmax(ind),model.dv(ind),1),mean(mm)) 30], [], [], [], [], x(:,1), x(:,2), nsolvent, 'silver', dagtop);
        p(:,j) = Simplex('CavityTransmissionFit', [polyval(polyfit(model.tmax(ind),model.dv(ind),1),mean(mm)) 30], [], [], [], [], x(:,1), x(:,2), nsolvent, [ag_lam ag_re+1i*ag_im], dagtop);
        num(j) = str2double(fnames(j).name(1:end-4));
        xlabel('wavelength (nm)'); ylabel('transmission (a.u.)');
        ax = axis;
        s = {['\itd\rm = ' mnum2str(p(1,j),3,1) ' nm'],['\itd\rm_{Ag} = ' mnum2str(p(2,j),3,1) ' nm']};
        text(ax(1)+0.7*diff(ax(1:2)),ax(3)+0.8*diff(ax(3:4)),s);
        eval(['print -dpng -r300 ' '''' dname fnames(j).name(1:end-4) '''' 'New'])
    end
    x = load([dname 'transmissions and decays.txt']);
    dv = zeros(size(x,1),2); dag = dv; tmax = zeros(size(x,1),1);
    for j=1:size(x,1)
        dv(j,:) = [p(1,num==x(j,1)) p(1,num==x(j,2))];
        dag(j,:) = [p(2,num==x(j,1)) p(2,num==x(j,2))];
        tmax(j)= x(j,3);
    end
	eval(['save ' '''' dname 'TransmissionFitNew' '''' ' tmax dv dag'])

    % modeling
    eval(['load ' '''' dname 'TransmissionFitNew' ''''])
    dv = mean(dv,2);
    dag = mean(dag,2);
    %model  = AlexeyModel('silver', nsolvent, lamex, dag, dagtop, dv);
    model  = AlexeyModel([ag_lam ag_re+1i*ag_im], nsolvent, lamex, dag, dagtop, dv);
    model.tmax = tmax;
    eval(['save ' '''' dname 'modelNew' '''' ' model'])
end
    
if 0
    ltdata = load('c:\Joerg\Doc\Nanophotonics&Optics\AlexeyChizhik\RotoDiffCavity\Glycerol R6G.txt');
    spectrum = load('c:\Joerg\Doc\Nanophotonics&Optics\AlexeyChizhik\RotoDiffCavity\spectrum R6G in glycerol.txt');
    load('c:\Joerg\Doc\Nanophotonics&Optics\AlexeyChizhik\RotoDiffCavity\Converted transmissions\R6G in glycerol 20111121\modelNew.mat','model');
    ind=ltdata(:,1)>=510 & ltdata(:,1)<=650;
    [parag taufitg lamg taug tau0g] = AlexeyCavityFitFull(ltdata(ind,:), spectrum, model, 1/1000);
    
    ltdata = load('c:\Joerg\Doc\Nanophotonics&Optics\AlexeyChizhik\RotoDiffCavity\Water R6G.txt');
    ltdata(19,:) = [];
    spectrum = load('c:\Joerg\Doc\Nanophotonics&Optics\AlexeyChizhik\RotoDiffCavity\spectrum R6G in water.txt');
    load('c:\Joerg\Doc\Nanophotonics&Optics\AlexeyChizhik\RotoDiffCavity\Converted transmissions\R6G in water 20111121\modelNew.mat','model');
    ind=ltdata(:,1)>=510 & ltdata(:,1)<=650;
    [paraw taufitw lamw tauw tau0w] = AlexeyCavityFitFull(ltdata(ind,:), spectrum, model, 1000);
    
    lam = linspace(min([lamw(1) lamg(1)]),max([lamw(end) lamg(end)]));
    plot(lamg,taufitg,lamw,taufitw,lamg,taug,'or',lamw,tauw,'ob',lam,tau0g*ones(size(lam)),':r',lam,tau0w*ones(size(lam)),':b')
    axis([505 655 1 4.5])
    xlabel('max. transmission wavelength (nm)')
    ylabel('lifetime (ns)')
    legend({'glycerol','water'},3)
end

if 0 % excitation figure for paper
   load('c:\Joerg\Doc\Nanophotonics&Optics\AlexeyChizhik\RotoDiffCavity\Converted transmissions\R6G in water 20111121\modelNew.mat')
   exc = model.exc(20);
   ind = exc.rho(:,1)<0.5;
   [~, ep] = FocusImage2D(exc.rho,exc.z,cat(4,cat(3,exc.fxc,exc.fxs),cat(3,exc.fyc,exc.fys)));
   ev = FocusImage2D(exc.rho,exc.z,cat(3,exc.fzc,exc.fzs));
   ev = ev/max(ep(:));
   ev = ev(ind,:);
   ep = ep/max(ep(:));
   ep = ep(ind,:);
   close all
   mpcolor([-flipud(exc.rho(ind,:));exc.rho(ind,:)]*1e3,1e3*[exc.z(ind,:);exc.z(ind,:)],[flipud(ep);ep])
   hold on
   mpcolor(exc.rho(ind,:)*1e3,1e3*exc.z(ind,:),ev)
   patch([-max(max(exc.rho(ind,:))) max(max(exc.rho(ind,:))) max(max(exc.rho(ind,:))) -max(max(exc.rho(ind,:)))]*1e3,1e3*exc.z(1) + 1e3*[0 0 -0.03 -0.03],[0.7 0.7 0.7],'edgecolor','none')
   patch([-max(max(exc.rho(ind,:))) max(max(exc.rho(ind,:))) max(max(exc.rho(ind,:))) -max(max(exc.rho(ind,:)))]*1e3,1e3*exc.z(end) + 1e3*[0 0 0.085 0.085],[0.7 0.7 0.7],'edgecolor','none')
   hold off
   axis image 
   %set(gca,'XAxisLocation','top')
   set(gca,'xticklabel',abs(get(gca,'xtick')))
   set(gca,'linewidth',0.1)
   h = colorbar;
   set(h,'linewidth',0.5)
   xlabel('\rho (nm)');
   ylabel('\itz \rm (nm)')
   print -dpng -r600 AlexeyExcitation
end

if 0 % emission figure for paper
   load('c:\Joerg\Doc\Nanophotonics&Optics\AlexeyChizhik\RotoDiffCavity\Converted transmissions\R6G in water 20111121\modelNew.mat')
   j = 20;
   lam = 500:650;
   zmax = 1e3*model.exc(j).z(end) + 0.5*model.exc(j).z(1);
   z0 = linspace(0,zmax,model.dz);
   z = linspace(0,zmax,100);
   lp = interp2(model.lambda,z0',model.lp(:,:,j),lam,z','spline');
   lv = interp2(model.lambda,z0',model.lv(:,:,j),lam,z','spline');
   close all
   tmp = contour(lam,[z max(z)+z],log10([lp;lv]),[0 0]);
   pcolor(lam,[z zmax+z],log10([lp;lv]))
   yt = get(gca,'ytick');
   yt = yt(yt<zmax);
   set(gca,'ytick',[yt zmax+yt])
   set(gca,'yticklabel',[yt yt])
   shading interp
   line(tmp(1,2:end),tmp(2,2:end),'color','w','linewidth',2)
   colormap jet
   line(lam,zmax*ones(size(lam)),'color','k','linewidth',1)
   set(gca,'linewidth',1)
   h = colorbar;
   set(h,'linewidth',1)
   set(gca,'DataAspectratio',[1 4 1])
   print -dpng -r600 AlexeyEmission
end

