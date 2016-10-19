% Confocal anisotropy factors
close all
clear all

if 1 % compute with S0->S1 optical saturation 
    NA = 1.2;
    fd = 3e3;
    n0 = 1.333;
    n = 1.333;
    n1 = 1.333;
    d0 = [];
    d = [];
    d1 = [];
    lamex = 0.64;
    overv = 1e3:250:5e3;
    focpos = 30;
    av = 150/2; % confocal pinhole radius
    lamem = 0.670;
    mag = 60;
    zpin = 0;
    atf = [];
    resolution = 20;
    kappa = 1;

    rhofield = [0 3]; % war [0 10]
    zfield = [25 35]; % war [20 40];

    gammav = (0:15:90)/180*pi;
    for jo=1:length(overv)
        over = overv(jo);
        exc = GaussExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, atf, resolution);
        ind = round(size(exc.z,2)/2);
        ff = abs(exc.fxc(:,ind,1)).^2 + abs(exc.fyc(:,ind,1)).^2 + abs(exc.fzc(:,ind,1)).^2;
        tmp = Simplex('Gauss',[0 0.1],[0 0],[0 inf],[],[],[-flipud(exc.rho(:,1)); exc.rho(:,1)],[flipud(ff); ff]');
        w0(jo) = 2*tmp(2);
        
        for jg = 1:length(gammav)
            gamma = gammav(jg);
            
            [aniso(jo,jg).ux aniso(jo,jg).uy aniso(jo,jg).pp aniso(jo,jg).po aniso(jo,jg).op aniso(jo,jg).oo] = GaussExc2RotoFCS(exc, lamem, mag, av, zpin, atf, gamma);
            
            eval(['save ConfoAniso12NA aniso w0 overv'])
        end
    end
end

if 0
    clear all
    load ConfoAniso12NA
    ww = w0;
    jmax = size(aniso,1);
    for j = 1:size(aniso,1)
        for k = 1:size(aniso,2)
            ux(j,k) = aniso(j,k).ux;
            uy(j,k) = aniso(j,k).uy;
            pp(:,:,j,k) = aniso(j,k).pp;
            po(:,:,j,k) = aniso(j,k).po;
            op(:,:,j,k) = aniso(j,k).op;
            oo(:,:,j,k) = aniso(j,k).oo;
        end
    end
    plot(0:90,interp1((0:15:90),((ux-uy)./(ux+2*uy))',0:90,'cubic'))
    for j = 1:size(aniso,1)
        for k = 1:size(aniso,2)
            op(:,:,j,k) = op(:,:,j,k)/pp(1,1,j,k);
            po(:,:,j,k) = po(:,:,j,k)/pp(1,1,j,k);
            oo(:,:,j,k) = oo(:,:,j,k)/pp(1,1,j,k);            
            pp(:,:,j,k) = pp(:,:,j,k)/pp(1,1,j,k);
        end
    end
    figure
    plot(0:90, interp1(0:15:90,squeeze(pp(1,1,:,:))',0:90,'cubic'), 'r', 0:90, interp1(0:15:90,squeeze(po(1,1,:,:))',0:90,'cubic'), 'g', ...
        0:90, interp1(0:15:90,squeeze(op(1,1,:,:))',0:90,'cubic'), 'b', 0:90, interp1(0:15:90,squeeze(oo(1,1,:,:))',0:90,'cubic'), 'c')

    over = overv;
    save ConfoAnisoCoefficients ux uy pp po op oo w0 over
end

if 0 % compute with S0->S1 optical saturation 
    NA = 1e-3;
    fd = 3e3;
    n0 = 1.333;
    n = 1.333;
    n1 = 1.333;
    d0 = [];
    d = [];
    d1 = [];
    lamex = 0.64;
    over = 1e3;
    focpos = 30;
    av = 150/2; % confocal pinhole radius
    lamem = 0.670;
    mag = 60;
    zpin = 0;
    atf = [];
    resolution = 20;
    kappa = 1;

    rhofield = [0 1];
    zfield = [29 31];

    exc = GaussExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, atf, resolution);
    [ux uy pp po op] = GaussExc2RotoFCS(exc, lamem, mag, av, zpin, atf);
    
    eval('save ConfoAniso0 ux uy pp po op')
end


if 0 % make Aniso figure Anastasia
    if ~(exist('aniso')==1)
        load ConfoAnisoSat12NA
    end
    jsat = 1;
    for jo=1:length(overv)
        for L=0:2:4
            tmpp(jo,L/2+1) = aniso(jsat,jo).pp(L+1,1) + 2*sum(aniso(jsat,jo).pp(L+1,2:end));
            tmpo(jo,L/2+1) = aniso(jsat,jo).po(L+1,1) + 2*sum(real(1i.^(1:6).*aniso(jsat,jo).po(L+1,2:end)));
            tmop(jo,L/2+1) = aniso(jsat,jo).op(L+1,1) + 2*sum(real(1i.^(1:6).*aniso(jsat,jo).op(L+1,2:end)));
        end
        tmpp(jo,:) = tmpp(jo,:)/aniso(jsat,jo).po(1,1);
        tmpo(jo,:) = tmpo(jo,:)/aniso(jsat,jo).po(1,1);
        tmop(jo,:) = tmop(jo,:)/aniso(jsat,jo).po(1,1);
    end
    x = overv/1e3;
    plot(x,3*tmpp(:,1),'r',x,3*tmpp(:,2),'r--',x,3*tmpp(:,3),'r:',x,tmpo(:,1),'b',x,tmpo(:,2),'b--',x,tmpo(:,3),'b:',x,9*tmop(:,1),'g',x,9*tmop(:,2),'g--',x,9*tmop(:,3),'g:')
    legend({'|| \times ||, \itl\rm = 0, \times3', '|| \times ||, \itl\rm = 2, \times3', '|| \times ||, \itl\rm = 4, \times3','|| \times \perp, \itl\rm = 0', '|| \times \perp, \itl\rm = 2', '|| \times \perp, \itl\rm = 4','\perp \times ||, \itl\rm = 0, \times9', '\perp \times ||, \itl\rm = 2, \times9', '\perp \times ||, \itl\rm = 4, \times9'},'Location','NorthEastOutside')    
    xlabel('laser beam diameter (mm)')
    ylabel('$$\Sigma_{m}\ u^{*}_{2,lm}\ u_{1,lm}$$','Interpreter','latex')
    grid
    
    figure
    drot = 1/20/6;
    t = 0:100;
    jo = length(overv);
    plot(t,tmpp(jo,1)+tmpp(jo,2)*exp(-6*drot*t)+tmpp(jo,3)*exp(-20*drot*t),...
        t,tmpo(jo,1)+tmpo(jo,2)*exp(-6*drot*t)+tmpo(jo,3)*exp(-20*drot*t),...
        t,tmop(jo,1)+tmop(jo,2)*exp(-6*drot*t)+tmop(jo,3)*exp(-20*drot*t))
    xlabel('time (ns)');
    ylabel('correlation (a.u.)')
    legend({'pp','po','op'})
end

if 0 % make Aniso figure Christoph
    if ~(exist('aniso')==1)
        load ConfoAnisoSat12NA
    end
    jsat = 1;
    for jo=1:length(overv)
        for L=0:2:4
            tmpp(jo,L/2+1) = aniso(jsat,jo).pp(L+1,1) + 2*sum(aniso(jsat,jo).pp(L+1,2:end));
            tmpo(jo,L/2+1) = aniso(jsat,jo).po(L+1,1) + 2*sum(real(aniso(jsat,jo).po(L+1,2:end)));
            tmop(jo,L/2+1) = aniso(jsat,jo).op(L+1,1) + 2*sum(real(aniso(jsat,jo).op(L+1,2:end)));
        end
        tmpp(jo,:) = tmpp(jo,:)/aniso(jsat,jo).po(1,1);
        tmpo(jo,:) = tmpo(jo,:)/aniso(jsat,jo).po(1,1);
        tmop(jo,:) = tmop(jo,:)/aniso(jsat,jo).po(1,1);
    end
    x = overv/1e3;
    plot(x,tmpo(:,1),'r',x,tmpo(:,2),'r--',x,tmpo(:,3),'r:',x,3*tmpp(:,1),'b',x,3*tmpp(:,2),'b--',x,3*tmpp(:,3),'b:',x,9*tmop(:,1),'g',x,9*tmop(:,2),'g--',x,9*tmop(:,3),'g:')
    legend({'|| \times ||, \itl\rm = 0', '|| \times ||, \itl\rm = 2', '|| \times ||, \itl\rm = 4', '|| \times \perp, \itl\rm = 0, \times3', '|| \times \perp, \itl\rm = 2, \times3', '|| \times \perp, \itl\rm = 4, \times3', '\perp \times \perp, \itl\rm = 0, \times9', '\perp \times \perp, \itl\rm = 2, \times9', '\perp \times \perp, \itl\rm = 4, \times9'},'Location','NorthEastOutside')    
    xlabel('laser beam diameter (mm)')
    ylabel('$$\Sigma_{m}\ u^{*}_{2,lm}\ u_{1,lm}$$','Interpreter','latex')
    grid
    
    figure
    drot = 1/20/6;
    t = 0:100;
    jo = length(overv);
    plot(t,tmpo(jo,1)+tmpo(jo,2)*exp(-6*drot*t)+tmpo(jo,3)*exp(-20*drot*t),...
        t,tmpp(jo,1)+tmpp(jo,2)*exp(-6*drot*t)+tmpp(jo,3)*exp(-20*drot*t),...
        t,tmop(jo,1)+tmop(jo,2)*exp(-6*drot*t)+tmop(jo,3)*exp(-20*drot*t))
    xlabel('time (ns)');
    ylabel('correlation (a.u.)')
    legend({'|| \times ||','|| \times \perp','\perp \times \perp'})
end

if 0 % make Aniso figure a la Rigler-Ehrenfest
    figure
    if ~(exist('aniso')==1)
        load ConfoAnisoSat12NA
    end
    jsat = 1;
    for jo=1:length(overv)
        for L=0:2:4
            tmpp_po(jo,L/2+1) = aniso(jsat,jo).pp(L+1,1) + 2*sum(aniso(jsat,jo).pp(L+1,2:end));
            tmpp_pp(jo,L/2+1) = aniso(jsat,jo).po(L+1,1) + 2*sum(aniso(jsat,jo).po(L+1,2:end));
            tmpo_po(jo,L/2+1) = aniso(jsat,jo).op(L+1,1) + 2*sum(aniso(jsat,jo).op(L+1,2:end));
            tmpp_op(jo,L/2+1) = aniso(jsat,jo).pp(L+1,1) + 2*sum(real(i.^(1:6).*aniso(jsat,jo).pp(L+1,2:end)));
            tmpo_op(jo,L/2+1) = aniso(jsat,jo).op(L+1,1) + 2*sum(real(i.^(1:6).*aniso(jsat,jo).op(L+1,2:end)));
            tmpp_oo(jo,L/2+1) = aniso(jsat,jo).po(L+1,1) + 2*sum(real(i.^(1:6).*aniso(jsat,jo).po(L+1,2:end)));
        end
        tmpp_pp(jo,:) = tmpp_pp(jo,:)/aniso(jsat,jo).po(1,1);
        tmpp_po(jo,:) = tmpp_po(jo,:)/aniso(jsat,jo).po(1,1);
        tmpo_po(jo,:) = tmpo_po(jo,:)/aniso(jsat,jo).po(1,1);
        tmpp_op(jo,:) = tmpp_op(jo,:)/aniso(jsat,jo).po(1,1);
        tmpo_op(jo,:) = tmpo_op(jo,:)/aniso(jsat,jo).po(1,1);
        tmpp_oo(jo,:) = tmpp_oo(jo,:)/aniso(jsat,jo).po(1,1);
    end
    x = overv/1e3;
    plot(x,tmpp_pp(:,1)+2*tmpp_po(:,1)+tmpo_po(:,1),'r',x,tmpp_pp(:,2)+2*tmpp_po(:,2)+tmpo_po(:,2),'r--',...
        x,2*tmpp_op(:,2)+tmpo_op(:,2)+tmpp_oo(:,2),'b--')
    legend({'|| \times || and || \times \perp, \itl\rm = 0', '|| \times ||, \itl\rm = 2', '|| \times \perp, \itl\rm = 2'},'Location','East')
    xlabel('laser beam diameter (mm)')
    ylabel('$$\Sigma_{m}\ u^{*}_{2,lm}\ u_{1,lm}$$','Interpreter','latex')
    grid

    figure
    drot = 1/20/6;
    t = 0:100;
    jo = length(overv);
    plot(t,tmpp_pp(jo,1)+2*tmpp_po(jo,1)+tmpo_po(jo,1) + (tmpp_pp(jo,2)+2*tmpp_po(jo,2)+tmpo_po(jo,2))*exp(-6*drot*t) + (tmpp_pp(jo,3)+2*tmpp_po(jo,3)+tmpo_po(jo,3))*exp(-20*drot*t),...
        t,2*tmpp_op(jo,1)+tmpo_op(jo,1)+tmpp_oo(jo,1) + (2*tmpp_op(jo,2)+tmpo_op(jo,2)+tmpp_oo(jo,2))*exp(-6*drot*t) + (2*tmpp_op(jo,3)+tmpo_op(jo,3)+tmpp_oo(jo,3))*exp(-20*drot*t))
    xlabel('time (ns)');
    ylabel('correlation (a.u.)')

end

if 0 % make ConfoAnisoSatMovie
    a=real(ux.*conj(uyr));
    b=real(ux.*conj(uy));
    for k=1:size(a,1) 
        a(k,1,:) = 0.5*a(k,1,:); 
        b(k,1,:) = 0.5*b(k,1,:); 
    end
    close; 
    for k=1:16; 
        mbar3(0:6,0:6,squeeze(b(:,:,k)./max(max(a(:,:,k)))./(abs(b(:,:,k))>0)),(0.5:7)'*ones(1,7)); 
        view([40 30]); 
        axis square; 
        axis tight; 
        axis([ax(1:4) -0.1 0.45]); 
        title(['\itI_{max}\rm = ' mnum2str(satv(k),1,1) ' \itI_{sat}']); 
        text(7.5,2.4,-0.1,'L'); 
        text(3.2,-1.8,-0.1,'M'); 
        eval(['print -dpng -r300 tmp' mint2str(k,2)]); 
    end
end

% t=0:0.1:100; tau1=40; tau2=6/20*tau1; 
% plot(t,2*(pp(1)+pp(2)*exp(-t/tau1)+pp(3)*exp(-t/tau2)),t,po(1)+po(2)*exp(
% -t/tau1)+po(3)*exp(-t/tau2),t,op(1)+op(2)*exp(-t/tau1)+op(3)*exp(-t/tau2))