function [tau, lam, data, irf] = AlexeyLifetimeRead(irfname,dyefolder)

    irf = load(irfname);
    irf(:,2) = abs(irf(:,2)-mean(irf(end-100:end,2)));

    if strcmp(dyefolder(end),'\')==1
        fnames = dir([dyefolder '*.txt']);
        for j=1:length(fnames)
            data(:,:,j) = load([dyefolder fnames(j).name]);
            lam(j) = str2num(fnames(j).name(1:3));
        end
        
        for j=1:size(data,3)
            [cx, tmp] = DistFluofit(irf(:,2), data(:,2,j), 12.5, 0.008, [-50 50]);
            tau(j) = cx'*tmp/sum(cx);
        end
        close all    
        plot(lam,tau,'o-');
        ylabel('mean excited state lifetime (ns)') 
        xlabel('tranmission max wavelength (nm)')
    else
       data = load(dyefolder);
       [cx, tau] = DistFluofit(irf(:,2), data(:,2), 12.5, 0.008, [-50 50]);
       tau = cx'*tau/sum(cx);
       lam = [];
    end
end

% tau0 = AlexeyLifetimeRead('c:\Joerg\Doc\Nanophotonics&Optics\AlexeyChizhik\RotoDiffCavity\Laser signal.txt','c:\Joerg\Doc\Nanophotonics&Optics\AlexeyChizhik\RotoDiffCavity\free space lifetime oregongreen.txt')