clear all
close all
dirname = 'c:\Joerg\Doc\Nanophotonics&Optics\AlexeyChizhik\';
dirname = 'c:\Joerg\Doc\Nanophotonics&Optics\AlexeyChizhik\Data2009_04_05\';

nglass = 1.518;
noil = 1; %1.518; 
nsio = 1.46;
npoly = 1.489;
dagbottomv = 30:60;
dagtopv = 2*dagbottomv;
dsiobottom = 30;
dpoly = 70;
doilv = 5:100;

theta = (0.5:200)/200*asin(1.3/nglass);

freespec = load('c:\Joerg\Doc\Nanophotonics&Optics\AlexeyChizhik\SMSpectra\Spectra2008Sep\pi freespace.txt','ascii');
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
                n0 = [nglass nagex nsio]; d0 = [dagbottom dsiobottom];
                n = npoly; d = dpoly;
                n1 = [noil nagex nglass]; d1 = [doil dagtop];
                [rp(jl,jd,ja), rs(jl,jd,ja), tp(jl,jd,ja), ts(jl,jd,ja)] = Fresnel(nglass,[n0 n n1],[d0 d d1]/lamem*2*pi);
            end
        end
    end
    save 'c:\Joerg\Doc\Nanophotonics&Optics\AlexeyChizhik\CavityTransmissionAir' lambda dagbottomv dagtopv doilv nglass noil nsio npoly dsiobottom dpoly rp rs tp ts
end

if 0 % Cavity lifetime modelling
    load metals
    for ja = 1:10:31;
        dagbottom = dagbottomv(ja);
        dagtop = dagtopv(ja);
        for jd = 1:length(doilv)
            doil = doilv(jd);
            for jw = 1:length(lambda)
                lamem = lambda(jw);
                nagem = interp1(wavelength,silver,lamem);
                n0 = [nglass nagem nsio]; d0 = [dagbottom dsiobottom];
                n = npoly; d = dpoly;
                n1 = [noil nagem nglass]; d1 = [doil dagtop];

                [lvd,lvu,lpd,lpu,qvd,qvu,qpd,qpu] = LifetimeL(0,n0,n,n1,d0/lamem*2*pi,d/lamem*2*pi,d1/lamem*2*pi);
                [v,pc,ps] = DipoleL(theta,0,n0,n,n1,d0/lamem*2*pi,d/lamem*2*pi,d1/lamem*2*pi);

                bv(jw,jd) = (sin(theta)*abs(v).^2)';
                bp(jw,jd) = (0.5*sin(theta)*(abs(pc).^2+abs(ps).^2))';
                lv(jw,jd) = 4/3*n./(qvd+qvu)';
                lp(jw,jd) = 4/3*n./(qpd+qpu)';
            end
            disp([ja jd])
        end
        save(['c:\Joerg\Doc\Nanophotonics&Optics\AlexeyChizhik\CavityLifetimeAir' mint2str(dagbottomv(ja),2) 'nm'],'lambda','dagbottom','dagtopv','doilv','nglass','noil','nsio','npoly','dsiobottom','dpoly','lv','lp','bv','bp')
    end
end


if 1 % lifetime fitting
    names = dir('c:\Joerg\Doc\Nanophotonics&Optics\AlexeyChizhik\Data2009_04_05\*.txt');
    load('c:\Joerg\Doc\Nanophotonics&Optics\AlexeyChizhik\CavityTransmissionAir','tp')
    
    ja = 21; % picks the 30 nm Ag layer
    [pos,pos] = max(abs(tp(:,:,ja)).^2);
    pos = lambda(pos);
    [pos, ind] = unique(pos);
    load(['c:\Joerg\Doc\Nanophotonics&Optics\AlexeyChizhik\CavityLifetimeAir' mint2str(dagbottomv(ja)) 'nm'])
    doilv = doilv(ind);
    radp = sum((emission'*ones(1,length(doilv)))./lp(:,ind))/sum(emission);
    radv = sum((emission'*ones(1,length(doilv)))./lv(:,ind))/sum(emission);
    
    for j=1:length(names)
        x = load([dirname names(j).name]);
        doil = interp1(pos,doilv,x(:,1));
        radi = interp1(doilv,radp,doil);
        c = polyfit(radi,(1./x(:,2)),1);
        tau0(j) = 1/sum(c);
        qy(j) = c(1)*tau0(j);
        plot(x(:,1),x(:,2),'o',x(:,1),1./polyval(c,radi))
        ax = axis;
        text(ax(1)+0.7*diff(ax(1:2)), ax(3)+0.7*diff(ax(3:4)), {['\tau_0 = ' mnum2str(tau0(j),1,1) ' ns'],['\phi_{\itf} = ' mnum2str(qy(j),1,1)]});
        xlabel('max. transmission wavelength (nm)');
        ylabel('lifetime (ns)');
        eval(['print -dpng -r300 ' '''' names(j).name(1:end-4) ' ' mint2str(dagbottomv(ja),2) 'nm'''])
    end
end


if 0 % slightly sophisticated lifetime fitting
    names = dir('c:\Joerg\Doc\Nanophotonics&Optics\AlexeyChizhik\Data2009_04_05\*.txt');
    load('c:\Joerg\Doc\Nanophotonics&Optics\AlexeyChizhik\CavityTransmissionAir','tp')
    
    ja = 21; % picks the 30 nm Ag layer
    [pos,pos] = max(abs(tp(:,:,ja)).^2);
    pos = lambda(pos);
    [pos, ind] = unique(pos);
    load(['c:\Joerg\Doc\Nanophotonics&Optics\AlexeyChizhik\CavityLifetimeAir' mint2str(dagbottomv(ja)) 'nm'])
    doilv = doilv(ind);
    radp = sum((emission'*ones(1,length(doilv)))./lp(:,ind))/sum(emission);
    radv = sum((emission'*ones(1,length(doilv)))./lv(:,ind))/sum(emission);
    
    for j=1:length(names)
        x = load([dirname names(j).name]);
        doil = interp1(pos,doilv,x(:,1));
        radi = interp1(doilv,radp,doil);
        p = Simplex('AffineFit',[1 1],[-20 1],[20 1],[],[], x(:,1), 1./x(:,2), pos, radp, [], 2, 0, 1);
        [err, c, z] = AffineFit(p, x(:,1), 1./x(:,2), pos, radp, [], 2, 0, 1);
        tau0(j) = 1/sum(c);
        qy(j) = c(2)*tau0(j);
        plot(x(:,1),x(:,2),'o',x(:,1),1./z)
        ax = axis;
        text(ax(1)+0.7*diff(ax(1:2)), ax(3)+0.7*diff(ax(3:4)), {['\tau_0 = ' mnum2str(tau0(j),1,1) ' ns'],['\phi_{\itf} = ' mnum2str(qy(j),1,1)]});
        xlabel('max. transmission wavelength (nm)');
        ylabel('lifetime (ns)');
        eval(['print -dpng -r300 ' '''' names(j).name(1:end-4) ' ' mint2str(dagbottomv(ja),2) 'nm' mint2str(round(p(1)),2) 'shift'''])
    end

end


if 0 % sophisticated lifetime fitting
    names = dir('c:\Joerg\Doc\Nanophotonics&Optics\AlexeyChizhik\Data2009_04_05\*.txt');
    load('c:\Joerg\Doc\Nanophotonics&Optics\AlexeyChizhik\CavityLifetimeAir50nm')
    
    radp = sum((emission'*ones(1,length(doilv)))./squeeze(lp))/sum(emission);
    %radv = sum((emission'*ones(1,length(doilv)))./squeeze(lv))/sum(emission);

    for j=1:length(names)
        x = load([dirname names(j).name]);
        close; 
        p = Simplex('AffineFit',[-466 1],[],[],[],[], x(:,1), 1./x(:,2), doilv, radp, [], 2, 0, 1);
        [err, c, z] = AffineFit(p, x(:,1), 1./x(:,2), doilv, radp, [], 2, 0, 1);
        tau0(j) = 1/sum(c);
        qy(j) = c(2)*tau0(j);
        plot(x(:,1),x(:,2),'o',x(:,1),1./z)
        ax = axis;
        text(ax(1)+0.7*diff(ax(1:2)), ax(3)+0.7*diff(ax(3:4)), {['\tau_0 = ' mnum2str(tau0(j),1,1) ' ns'],['\phi_{\itf} = ' mnum2str(qy(j),1,1)]});
        xlabel('max. transmission wavelength (nm)');
        ylabel('lifetime (ns)');
        eval(['print -dpng -r300 ' '''' names(j).name(1:end-4) ' fit50nm'''])
    end
end
