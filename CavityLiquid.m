clear all
close all

NA = 1.25;
nglass = 1.518;
nsolvent = 1.5; % toluene
% nsolvent = 1.45; % ethylene glycol
dagbottomv = 30:80;
dsolventv = 20:1:140;

theta = (0.5:200)/200*asin(NA/nglass);

load('c:\Joerg\Doc\Nanophotonics&Optics\FrankSchleifenbaum\Liquid Cavity\PI_in_liquid_cavity.mat');
lambda = 450:0.5:750;
transmission = absorption;

if 0 % transmission spectra calculation
    load metals
    tp = zeros(length(lambda),length(dsolventv),length(dagbottomv));
    rp = tp;
    ts = tp;
    rs = tp;
    for ja = 1:length(dagbottomv)
        dagbottom = dagbottomv(ja);
        daubottom = 4/30*dagbottom;
        dautop = daubottom;
        dagtop = 2*dagbottom;
        for jd = 1:length(dsolventv)
            dsolvent = dsolventv(jd);
            for jl = 1:length(lambda)
                lamem = lambda(jl);
                nagex = interp1(wavelength,silver,lamem);
                nauex = interp1(wavelength,gold,lamem);
                n0 = [nglass nagex nauex]; d0 = [dagbottom daubottom ];
                n = nsolvent; d = dsolvent;
                n1 = [nauex nagex nglass]; d1 = [dautop dagtop];
                [rp(jl,jd,ja), rs(jl,jd,ja), tp(jl,jd,ja), ts(jl,jd,ja)] = Fresnel(nglass,[n0 n n1],[d0 d d1]/lamem*2*pi);
            end
        end
    end
    save TolueneCavityTransmission lambda dagbottomv nglass nsolvent rp rs tp ts
end

if 0 % transmission spectra fitting
    dirname = 'c:\Joerg\Doc\Nanophotonics&Optics\FrankSchleifenbaum\Liquid Cavity\';
    %load([dirname 'TolueneCavityTransmission'])
    [ind,ind] = max(transmission(:,2));
    ind = abs(transmission(:,1)-transmission(ind,1))<=15;
    errtot = zeros(1,length(dagbottomv));
    for ja=1:length(dagbottomv)
        for jd=1:length(dsolventv)
            c = [ones(sum(ind),1) cumsum(ones(sum(ind),1)) interp1(lambda,abs(tp(:,jd,ja)).^2,transmission(ind,1),'cubic')]\transmission(ind,2);
            tmp = [ones(sum(ind),1) cumsum(ones(sum(ind),1)) interp1(lambda,abs(tp(:,jd,ja)).^2,transmission(ind,1),'cubic')]*c;
            err(jd) = sum((transmission(ind,2)-tmp).^2.)/sum(transmission(ind,2));
        end
        [jd,jd] = min(err);
        c = [ones(sum(ind),1) cumsum(ones(sum(ind),1)) interp1(lambda,abs(tp(:,jd,ja)).^2,transmission(ind,1),'cubic')]\transmission(ind,2);
        modspec{ja} = [ones(sum(ind),1) cumsum(ones(sum(ind),1)) interp1(lambda,abs(tp(:,jd,ja)).^2,transmission(ind,1),'cubic')]*c;
        lamspec{ja} = transmission(ind,1);
        jdspec(ja) = jd;
        errtot(ja) = min(err);
    end
    [ja,ja] = min(errtot);
    plot(transmission(:,1), transmission(:,2)./max(transmission(:,2)),'o');
    hold on
    plot(lamspec{ja}, modspec{ja}/max(transmission(:,2)), 'm');
    hold off
    legend(num2str(dsolventv(jdspec(ja))'))
    xlabel('wavelength (nm)'); 
    ylabel('transmissivity (a.u.)'); 
    %print -dpng -r300 TransmissionSpec
    save TolueneCavityTransmissionFit ja jdspec modspec lamspec transmission
    disp([dagbottomv(ja) dsolventv(jdspec(ja))])
end

if 1 % Cavity lifetime modelling & excitation intensity toluene
    ja = 1; 
    dagbottom = dagbottomv(ja);
    dagtop = 2*dagbottom;
    daubottom = 4/30*dagbottom;
    dautop = daubottom;
    %load('c:\Joerg\Doc\Nanophotonics&Optics\FrankSchleifenbaum\Liquid Cavity\TolueneCavityTransmission.mat')
    %[ind,ind]= max(abs(tp(:,:,ja)));
    %c = polyfit(lambda(ind(lambda(ind)>450)),dsolventv(lambda(ind)>450),2);
    %dsolventv = polyval(c,lifetime1(:,1)*2*nsolvent+50);
    
    lambda = 450:10:800;
    dsolventv = 60:5:300;
    
    load metals
    % emission
    for jd = 1:length(dsolventv)
        dsolvent = dsolventv(jd);
        for jw = 1:length(lambda)
            lamem = lambda(jw);
            nagem = interp1(wavelength,silver,lamem);
            nauem = interp1(wavelength,gold,lamem);
            n0 = [nglass nagem nauem]; d0 = [dagbottom daubottom ];
            n = nsolvent; d = dsolvent;
            n1 = [nauem nagem nglass]; d1 = [dautop dagtop];
            zv = d/50:d/25:d;
            
            [lvd,lvu,lpd,lpu,qvd,qvu,qpd,qpu] = LifetimeL(zv/lamem*2*pi,n0,n,n1,d0/lamem*2*pi,d/lamem*2*pi,d1/lamem*2*pi);
            [v,pc,ps] = DipoleL(theta,zv/lamem*2*pi,n0,n,n1,d0/lamem*2*pi,d/lamem*2*pi,d1/lamem*2*pi);
            
            bv(:,jw,jd) = (sin(theta)*abs(v).^2)';
            bp(:,jw,jd) = (0.5*sin(theta)*(abs(pc).^2+abs(ps).^2))';
            lv(:,jw,jd) = 4/3*n./(qvd+qvu)';
            lp(:,jw,jd) = 4/3*n./(qpd+qpu)';
        end
        disp(jd)
    end

    % excitation
    lamex = 470;
    nagex = interp1(wavelength,silver,lamex);
    nauex = interp1(wavelength,gold,lamex);
    rhofield = [0 1];
    fd = 164.5e3/100;
    over = 3e3;
    focpos = 0;
    atf = [];
    n0 = [nglass nagex nauex]; d0 = [dagbottom daubottom]/1e3;
    n = nsolvent;
    n1 = [nauex nagex nglass]; d1 = [dautop dagtop]/1e3;
    for jd = 1:length(dsolventv)
        d = dsolventv(jd)/1e3;
        resolution = [20 lamex/1e3/(d/25)];
        exc(jd) = GaussExc(rhofield, [0, d], NA, fd, n0, n, n1, d0, d, d1, lamex/1e3, over, focpos, atf, resolution);
        disp(jd)
    end
    
    save('c:\Joerg\Doc\Nanophotonics&Optics\FrankSchleifenbaum\Liquid Cavity\TolueneCavityFluorescence030.mat', ...
        'lambda', 'zv', 'dsolventv', 'lv', 'lp', 'bv', 'bp', 'exc')
end

if 1 % Cavity lifetime modelling & excitation intensity glycol
    nsolvent = 1.45;
    ja = 1; 
    dagbottom = dagbottomv(ja);
    dagtop = 2*dagbottom;
    daubottom = 4/30*dagbottom;
    dautop = daubottom;
    %load('c:\Joerg\Doc\Nanophotonics&Optics\FrankSchleifenbaum\Liquid Cavity\TolueneCavityTransmission.mat')
    %[ind,ind]= max(abs(tp(:,:,ja)));
    %c = polyfit(lambda(ind(lambda(ind)>450)),dsolventv(lambda(ind)>450),2);
    %dsolventv = polyval(c,lifetime1(:,1)*2*nsolvent+50);
    
    lambda = 450:10:800;
    dsolventv = 60:5:300;
    
    load metals
    % emission
    for jd = 1:length(dsolventv)
        dsolvent = dsolventv(jd);
        for jw = 1:length(lambda)
            lamem = lambda(jw);
            nagem = interp1(wavelength,silver,lamem);
            nauem = interp1(wavelength,gold,lamem);
            n0 = [nglass nagem nauem]; d0 = [dagbottom daubottom ];
            n = nsolvent; d = dsolvent;
            n1 = [nauem nagem nglass]; d1 = [dautop dagtop];
            zv = d/50:d/25:d;
            
            [lvd,lvu,lpd,lpu,qvd,qvu,qpd,qpu] = LifetimeL(zv/lamem*2*pi,n0,n,n1,d0/lamem*2*pi,d/lamem*2*pi,d1/lamem*2*pi);
            [v,pc,ps] = DipoleL(theta,zv/lamem*2*pi,n0,n,n1,d0/lamem*2*pi,d/lamem*2*pi,d1/lamem*2*pi);
            
            bv(:,jw,jd) = (sin(theta)*abs(v).^2)';
            bp(:,jw,jd) = (0.5*sin(theta)*(abs(pc).^2+abs(ps).^2))';
            lv(:,jw,jd) = 4/3*n./(qvd+qvu)';
            lp(:,jw,jd) = 4/3*n./(qpd+qpu)';
        end
        disp(jd)
    end

    % excitation
    lamex = 470;
    nagex = interp1(wavelength,silver,lamex);
    nauex = interp1(wavelength,gold,lamex);
    rhofield = [0 1];
    fd = 164.5e3/100;
    over = 3e3;
    focpos = 0;
    atf = [];
    n0 = [nglass nagex nauex]; d0 = [dagbottom daubottom]/1e3;
    n = nsolvent;
    n1 = [nauex nagex nglass]; d1 = [dautop dagtop]/1e3;
    for jd = 1:length(dsolventv)
        d = dsolventv(jd)/1e3;
        resolution = [20 lamex/1e3/(d/25)];
        exc(jd) = GaussExc(rhofield, [0, d], NA, fd, n0, n, n1, d0, d, d1, lamex/1e3, over, focpos, atf, resolution);
        disp(jd)
    end
    
    save('c:\Joerg\Doc\Nanophotonics&Optics\FrankSchleifenbaum\Liquid Cavity\GlycolCavityFluorescence030.mat', ...
        'lambda', 'zv', 'dsolventv', 'lv', 'lp', 'bv', 'bp', 'exc')
end

if 1 % lifetime fitting
    load('c:\Joerg\Doc\Nanophotonics&Optics\FrankSchleifenbaum\Liquid Cavity\PI_in_liquid_cavity.mat')
    load('c:\Joerg\Doc\Nanophotonics&Optics\FrankSchleifenbaum\Liquid Cavity\TolueneCavityFluorescence030.mat')
    ind = emission(1,3)<lambda & lambda<emission(end,3);
    lam = lambda(ind);
    spec = interp1(emission(:,3),abs(emission(:,4)-mean(emission(1:10,4))),lam,'cubic'); % toluene
    spec = interp1(emission(:,1),abs(emission(:,2)-mean(emission(1:10,2))),lam,'cubic'); % toluene
    lifetime = lifetime1; % toluene

    for k=1:length(dsolventv)
        fp = abs(exc(k).fxc(:,:,1)).^2+0.5*sum(abs(exc(k).fxc(:,:,2:end)).^2,3)+0.5*sum(abs(exc(k).fxs).^2,3) + ...
            abs(exc(k).fyc(:,:,1)).^2+0.5*sum(abs(exc(k).fyc(:,:,2:end)).^2,3)+0.5*sum(abs(exc(k).fys).^2,3);
        fv = abs(exc(k).fzc(:,:,1)).^2+0.5*sum(abs(exc(k).fzc(:,:,2:end)).^2,3)+0.5*sum(abs(exc(k).fzs).^2,3);
        fp = exc(k).rho(:,1)'*fp/2;
        fv = exc(k).rho(:,1)'*fv;        
        feld = fp+fv;
        
        dp = spec*interp1(lambda, bp(:,:,k)', lam, 'cubic');
        dv = spec*interp1(lambda, bv(:,:,k)', lam, 'cubic');
        sp = spec*interp1(lambda, 1./lp(:,:,k)', lam, 'cubic');
        sv = spec*interp1(lambda, 1./lv(:,:,k)', lam, 'cubic');
        weight = real((2*((dp - dv).*sqrt(sp).*sqrt(sp - sv) + (dv.*sp - dp.*sv).*atanh(sqrt(sp - sv)./sqrt(sp))))...
            ./(sqrt(sp).*(sp - sv).^(3/2))).*feld;
        stot = (2*sp+sv)/sum(spec)/3;

        tauf(k) = weight*(1./stot)'/sum(weight(:));


        weight = real(2.*(sqrt(sp).*sqrt(sp - sv).*(dp.*(2*fp.*sp + fv.*sp - 5*fp.*sv + 2*fv.*sv) + ...
            dv.*(fv.*(-4*sp + sv) + fp.*(sp + 2*sv))) + 3*(dv.*sp - dp.*sv).*(fv.*sp - fp.*sv).*atanh(sqrt(sp - sv)./sqrt(sp)))./ ...
            (3*sqrt(sp).*(sp - sv).^(5/2)));
        
        stot = real((sqrt(sp).*sqrt(sp - sv).*(dp.*sv.*(-3*fv.*sp + fp.*(2*sp + sv)) + ...
            dv.*sp.*(-3*fp.*sv + fv.*(sp + 2*sv))) + sv.*(dv.*sp.*(-3*fv.*sp + fp.*(2*sp + sv)) + ...
            dp.*(fp.*sv.*(-4*sp + sv) + fv.*sp.*(2*sp + sv))).*atanh(sqrt(sp - sv)./sqrt(sp)))./(sp.^(3/2).*(sp - sv).^(5/2).*sv));
 
        taus(k) = (sum(spec)./weight)*stot';
    end
    close; 
    p1 = Simplex('CavityFit',[1 0 0.8],[-inf -inf 0],[inf inf 1],[],[],dsolventv,[tauf' taus'],lifetime(:,1),lifetime(:,2));
    [err,c1] = CavityFit(p1,dsolventv,[tauf' taus'],lifetime(:,1),lifetime(:,2));
    tauf1 = tauf; taus1 = taus;


    load('c:\Joerg\Doc\Nanophotonics&Optics\FrankSchleifenbaum\Liquid Cavity\PI_in_liquid_cavity.mat')
    load('c:\Joerg\Doc\Nanophotonics&Optics\FrankSchleifenbaum\Liquid Cavity\GlycolCavityFluorescence030.mat')
    ind = emission(1,3)<lambda & lambda<emission(end,3);
    lam = lambda(ind);
    spec = interp1(emission(:,1),abs(emission(:,2)-mean(emission(1:10,2))),lam,'cubic'); % glycol
    lifetime = lifetime2; % glycol

    for k=1:length(dsolventv)
        feld = abs(exc(k).fxc(:,:,1)).^2+0.5*sum(abs(exc(k).fxc(:,:,2:end)).^2,3)+0.5*sum(abs(exc(k).fxs).^2,3) + ...
            abs(exc(k).fyc(:,:,1)).^2+0.5*sum(abs(exc(k).fyc(:,:,2:end)).^2,3)+0.5*sum(abs(exc(k).fys).^2,3)+ ...
            abs(exc(k).fzc(:,:,1)).^2+0.5*sum(abs(exc(k).fzc(:,:,2:end)).^2,3)+0.5*sum(abs(exc(k).fzs).^2,3);
        feld = exc(k).rho(:,1)'*feld;
        
        dp = spec*interp1(lambda, bp(:,:,k)', lam, 'cubic');
        dv = spec*interp1(lambda, bv(:,:,k)', lam, 'cubic');
        sp = spec*interp1(lambda, 1./lp(:,:,k)', lam, 'cubic');
        sv = spec*interp1(lambda, 1./lv(:,:,k)', lam, 'cubic');
        cdet = real((2*((dp - dv).*sqrt(sp).*sqrt(sp - sv) + (dv.*sp - dp.*sv).*atanh(sqrt(sp - sv)./sqrt(sp))))...
            ./(sqrt(sp).*(sp - sv).^(3/2)));
        
        stot = (2*sp+sv)/sum(spec)/3;

        weight = cdet.*feld;
        weight = weight/sum(weight(:));
        tauf(k) = weight*(1./stot)';

        
        fp = abs(exc(k).fxc(:,:,1)).^2+0.5*sum(abs(exc(k).fxc(:,:,2:end)).^2,3)+0.5*sum(abs(exc(k).fxs).^2,3) + ...
            abs(exc(k).fyc(:,:,1)).^2+0.5*sum(abs(exc(k).fyc(:,:,2:end)).^2,3)+0.5*sum(abs(exc(k).fys).^2,3);
        fv = abs(exc(k).fzc(:,:,1)).^2+0.5*sum(abs(exc(k).fzc(:,:,2:end)).^2,3)+0.5*sum(abs(exc(k).fzs).^2,3);
        fp = exc(k).rho(:,1)'*fp;
        fv = exc(k).rho(:,1)'*fv;        
        
        dp = (spec'*ones(1,length(zv))).*interp1(lambda, bp(:,:,k)', lam, 'cubic').*(ones(length(spec),1)*fp);
        dv = (spec'*ones(1,length(zv))).*interp1(lambda, bv(:,:,k)', lam, 'cubic').*(ones(length(spec),1)*fv);        
        sp = interp1(lambda, 1./lp(:,:,k)', lam, 'cubic');
        sv = interp1(lambda, 1./lv(:,:,k)', lam, 'cubic');

        tmp = real((2*((dp - dv).*sqrt(sp).*sqrt(sp - sv) + (dv.*sp - dp.*sv).*atanh(sqrt(sp - sv)./sqrt(sp))))...
            ./(sqrt(sp).*(sp - sv).^(3/2)));
        
        dp = dp.*interp1(lambda,2*lp(:,:,k)',lam,'cubic');
        dv = dv.*interp1(lambda,lv(:,:,k)',lam,'cubic');        

        taus(k) = sum(sum(real((2*((dp - dv).*sqrt(sp).*sqrt(sp - sv) + (dv.*sp - dp.*sv).*atanh(sqrt(sp - sv)./sqrt(sp))))...
            ./(sqrt(sp).*(sp - sv).^(3/2)))))/sum(sum(tmp))/3;
    end
    close; 
    p2 = Simplex('CavityFit',[1 0 0.8],[-inf -inf 0],[inf inf 1],[],[],dsolventv,[tauf' taus'],lifetime(:,1),lifetime(:,2));
    [err,c2] = CavityFit(p2,dsolventv,[tauf' taus'],lifetime(:,1),lifetime(:,2));

    plot(lifetime1(:,1),lifetime1(:,2),'o',lifetime2(:,1),lifetime2(:,2),'o',...
        polyval(p1(1:end-1),dsolventv),c1(1)./(p1(end)./tauf1+1-p1(end))+c1(2)./(p1(end)./taus1+1-p1(end)),...
        polyval(p2(1:end-1),dsolventv),c2(1)./(p2(end)./tauf+1-p2(end))+c2(2)./(p2(end)./taus+1-p2(end)));

    plot(interp1(polyval(p1(1:end-1),dsolventv),dsolventv,lifetime1(:,1),'cubic'),lifetime1(:,2),'o',...
        interp1(polyval(p2(1:end-1),dsolventv),dsolventv,lifetime2(:,1),'cubic'),lifetime2(:,2),'o',...
        dsolventv,c1(1)./(p1(3)./tauf1+1-p1(3))+c1(2)./(p1(3)./taus1+1-p1(3)),...
        dsolventv,c2(1)./(p2(3)./tauf+1-p2(3))+c2(2)./(p2(3)./taus+1-p2(3)));
end
