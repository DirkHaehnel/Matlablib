clear all
close all

nglass = 1.518;
noil = 1.518; 
nsio = 1.46;
npoly = 1.489;
dagbottomv = 30:60;
dagtopv = 2*dagbottomv;
dsiobottom = 30;
dsiotop = 8;
dpoly = 70;
doilv = 10:0.1:60;
zv = 0:10:50; % dipole position within polymer layer

theta = (0.5:200)/200*asin(1.3/nglass);

% load c:\Joerg\Doc\Nanophotonics&Optics\FrankSchleifenbaum\Alexei\PISpectrum.mat
% lambda = mean(reshape(lambda,8,2^7));
% emission = mean(reshape(emission,8,2^7));
freespec = load('c:\Joerg\Doc\Nanophotonics&Optics\FrankSchleifenbaum\Alexei\Spectra2008Sep\pi freespace.txt','ascii');
lambda = mean(reshape(freespec(:,1),5,2*134));
emission = mean(reshape(freespec(:,2),5,2*134));

if 0 % transmission spectra calculation
    load metals
    tp = zeros(length(lambda),length(doilv),length(dagbottomv));
    rp = tp;
    ts = tp;
    rs = tp;
    for ja = 1:length(dagbottomv)
        dagbottom = dagbottomv(ja);
        dagtop = dagtopv(ja);
        for jd = 1:length(doilv)
            doil = doilv(jd);
            for jl = 1:length(lambda)
                lamem = lambda(jl);
                nagex = interp1(wavelength,silver,lamem);
                n0 = [nglass nagex nsio]; d0 = [dagbottom dsiobottom ];
                n = npoly; d = dpoly;
                n1 = [nsio noil nagex nglass]; d1 = [dsiotop doil dagtop];
                [rp(jl,jd,ja), rs(jl,jd,ja), tp(jl,jd,ja), ts(jl,jd,ja)] = Fresnel(nglass,[n0 n n1],[d0 d d1]/lamem*2*pi);
            end
        end
    end
    save CavityTransmission lambda dagbottomv dagtopv doilv nglass noil nsio npoly dsiobottom dsiotop dpoly rp rs tp ts
end

if 0 % transmission spectra fitting
    dirname = 'c:\Joerg\Doc\Nanophotonics&Optics\FrankSchleifenbaum\Alexei\Spectra2008Sep\';
    names = dir([dirname 'trans*']);
    load CavityTransmission
    errtot = zeros(1,length(dagbottomv));
    for ja=1:length(dagbottomv)
        for k=1:length(names)
            transspec = load([dirname names(k).name],'ascii');
            [ind,ind] = max(transspec(:,2));
            ind = abs(transspec(:,1)-transspec(ind,1))<=25;
            for jd=1:length(doilv)
                c = [ones(sum(ind),1) interp1(lambda,abs(tp(:,jd,ja)).^2,transspec(ind,1),'cubic')]\transspec(ind,2);
                tmp = [ones(sum(ind),1) interp1(lambda,abs(tp(:,jd,ja)).^2,transspec(ind,1),'cubic')]*c;
                err(jd) = sum((transspec(ind,2)-tmp).^2.)/sum(transspec(ind,2));
            end
            [jd,jd] = min(err);
            c = [ones(sum(ind),1) interp1(lambda,abs(tp(:,jd,ja)).^2,transspec(ind,1),'cubic')]\transspec(ind,2);
            modspec{k,ja} = [ones(sum(ind),1) interp1(lambda,abs(tp(:,jd,ja)).^2,transspec(ind,1),'cubic')]*c;
            lamspec{k,ja} = transspec(ind,1);
            jdspec(k,ja) = jd;
            errtot(ja) = errtot(ja) + min(err);
        end
    end
    [ja,ja] = min(errtot);
    k = 1;
    transspec = load([dirname names(k).name],'ascii');
    for k=2:length(names)
        tmp = load([dirname names(k).name],'ascii');
        transspec(:,end+1) = tmp(:,2);
    end
    plot(transspec(:,1), transspec(:,2:end)./(ones(size(transspec,1),1)*max(transspec(:,2:end))),'o');
    hold on
    for k=1:length(names)
        plot(lamspec{k,ja}, modspec{k,ja}/max(transspec(:,k+1)), 'm');
    end
    hold off
    legend(num2str(doilv(jdspec(:,ja))'))
    xlabel('wavelength (nm)'); 
    ylabel('transmissivity (a.u.)'); 
    %print -dpng -r300 TransmissionSpec
    
    save CavityTransmissionFit ja jdspec modspec lamspec transspec
end

if 0 % Cavity lifetime modelling
    load CavityTransmissionFit
    dagbottom = dagbottomv(ja);
    dagtop = 2*dagbottom;
    %doilv = doilv(jdspec(:,ja));
    doilv = 10:0.5:60;
    
    load metals
    for jd = 1:length(doilv)
        doil = doilv(jd);
        for jw = 1:length(lambda)
            lamem = lambda(jw);
            nagem = interp1(wavelength,silver,lamem);
            n0 = [nglass nagem nsio]; d0 = [dagbottom dsiobottom];
            n = npoly; d = dpoly;
            n1 = [nsio noil nagem nglass]; d1 = [dsiotop doil dagtop];
                
            [lvd,lvu,lpd,lpu,qvd,qvu,qpd,qpu] = LifetimeL(zv/lamem*2*pi,n0,n,n1,d0/lamem*2*pi,d/lamem*2*pi,d1/lamem*2*pi);
            [v,pc,ps] = DipoleL(theta,zv/lamem*2*pi,n0,n,n1,d0/lamem*2*pi,d/lamem*2*pi,d1/lamem*2*pi);
            
            bv(:,jw,jd) = (sin(theta)*abs(v).^2)';
            bp(:,jw,jd) = (0.5*sin(theta)*(abs(pc).^2+abs(ps).^2))';
            lv(:,jw,jd) = 4/3*n./(qvd+qvu)';
            lp(:,jw,jd) = 4/3*n./(qpd+qpu)';
        end
        disp(jd)
    end
    save CavityFluorescenceAll lambda emission zv doilv lv lp bv bp
end

if 1 % spectra fitting
    dirname = 'c:\Joerg\Doc\Nanophotonics&Optics\FrankSchleifenbaum\Alexei\Spectra2008Sep\';
    names = dir([dirname 'spec*']);
    load CavityFluorescenceAll
    load CavityTransmissionFit

    for k=1:length(names)
        fluospec = load([dirname names(k).name],'ascii');
        fluospec = fluospec(fluospec(:,1)>lambda(1) & fluospec(:,1)<lambda(end),:);
        
        for j=1:size(bp,3)
            para = squeeze(bp(:,:,j)./lp(:,:,j)).*(ones(size(bp,1),1)*emission);
            vert = squeeze(bv(:,:,j)./lv(:,:,j)).*(ones(size(bp,1),1)*emission);

            for jj=1:15
                tmp = interp1(lambda-jj+1,para',fluospec(:,1),'cubic',0);
                tmv = interp1(lambda-jj+1,vert',fluospec(:,1),'cubic',0);
                tmp = tmp./(ones(size(tmp,1),1)*sum(tmv));
                tmv = tmv./(ones(size(tmv,1),1)*sum(tmv));

                c = lsqnonneg([ones(size(fluospec,1),1) tmp tmv],fluospec(:,2) - min(fluospec(:,2)))
                tmp = [ones(size(fluospec,1),1) tmp tmv]*c+min(fluospec(:,2));
                err(j,jj) = sum((fluospec(:,2)-tmp).^2);
            end
        end
        [tmp,j] = min(err);
        [jj,jj] = min(tmp);
        j = j(jj);
        fac = 1+(jj-1)/500;
        para = squeeze(bp(:,:,j)./lp(:,:,j)).*(ones(size(bp,1),1)*emission);
        vert = squeeze(bv(:,:,j)./lv(:,:,j)).*(ones(size(bp,1),1)*emission);

        tmp = interp1(lambda-jj+1,para',fluospec(:,1),'cubic',0);
        tmv = interp1(lambda-jj+1,vert',fluospec(:,1),'cubic',0);
        tmp = tmp./(ones(size(tmp,1),1)*sum(tmv));
        tmv = tmv./(ones(size(tmv,1),1)*sum(tmv));

        c = lsqnonneg([ones(size(fluospec,1),1) tmp tmv],fluospec(1:end,2)-min(fluospec(:,2)))
        dat{k} = fluospec(1:end,2);
        fit{k} = [ones(size(fluospec,1),1) tmp tmv]*c+min(fluospec(:,2));
        jres(:,k) = [j;jj];
        cres(:,k) = c;

        area(lambda,emission/max(emission),'FaceColor',[0.8 0.8 0.8],'EdgeColor','none');
        hold on
        mm = max(fit{k}-mean(fit{k}(end-30:end)));
%         plot(jj+fluospec(:,1)-1,(fluospec(:,2)-mean(fit{k}(end-30:end)))/mm,'o',...
%             jj+fluospec(:,1)-1,(fit{k}-mean(fit{k}(end-30:end)))/mm);
        plot(fluospec(:,1)-1,(fluospec(:,2)-mean(fit{k}(end-30:end)))/mm,'o',...
            fluospec(:,1)-1,(fit{k}-mean(fit{k}(end-30:end)))/mm);
        hold off
        axis([500 750 -0.2 1.2])
        xlabel('wavelength (nm)');
        ylabel('emission intensity (a.u.)');
        text(630,1,['oil layer thickness = ' mnum2str(doilv(jres(1,k)),2,1) ' nm']);
        text(630,0.8,['spectral shift = ' int2str(-jj) ' nm']);
        %eval(['print -dpng -r300 ''PI-SM-Spec''' mint2str(k,1)'])
        drawnow
        covcor(k) = sum(fluospec(:,2).*fit{1})/sqrt(sum(fluospec(:,2).^2))/sqrt(sum(fit{k}.^2));
        %pause
    end
    cres(2:end,:)./(ones(12,1)*sum(cres(2:end,:)))
end

if 0 % lifetime
    for j=1:length(zv)
        tmp(j,:) = sum((emission'*ones(1,length(doilv)))./squeeze(lp(j,:,:)))/sum(emission);
        tmv(j,:) = sum((emission'*ones(1,length(doilv)))./squeeze(lv(j,:,:)))/sum(emission); 
    end
    figure(gcf); 
    for z=1:max(zv) 
        plot(doilv, interp1(zv,tmp,z,'cubic'), doilv, interp1(zv,tmv,z,'cubic')); 
        axis([10 150 0 10]); 
        text(60,9,['\itz_{dipole}\rm = ' mint2str(z,2) ' nm']); 
        xlabel('oil layer thickness (nm)'); 
        ylabel('rel. emission rate'); 
        legend({'horizontal','vertical'}); 
        drawnow
        eval(['print -dpng -r300 tmp' mint2str(z,2)])
    end    
end

if 0 % brightness
    for j=1:length(zv)
        tmp(j,:) = sum((emission'*ones(1,length(doilv))).*squeeze(bp(j,:,:)./lp(j,:,:)));
        tmv(j,:) = sum((emission'*ones(1,length(doilv))).*squeeze(bv(j,:,:)./lv(j,:,:))); 
    end
    figure(gcf); 
    for z=1:max(zv) 
        plot(doilv, interp1(zv,tmp,z,'cubic'), doilv, interp1(zv,tmv,z,'cubic')); 
        axis([10 150 0 4e5]); 
        text(60,3.8e5,['\itz_{dipole}\rm = ' mint2str(z,2) ' nm']); 
        xlabel('oil layer thickness (nm)'); 
        ylabel('emission intensity (a.u.)'); 
        legend({'horizontal','vertical'}); 
        drawnow
        eval(['print -dpng -r300 tmp' mint2str(z,2)])
    end    
end