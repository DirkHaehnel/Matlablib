clear all
close all

r2 = 40;
t = 0:1e3;
prn = 1;

if 0 % influence of depolarization/focusing
    load ConfoAnisoCoefficients w0
    
    wv = fliplr(linspace(min(w0),max(w0),300));
    
    z = zeros(length(t),4,length(wv));
    for jw=1:length(wv)
        z(:,:,jw) = RotoDiffCWFit(r2, t, [], [], [], [], wv(jw));
        p(jw) = Simplex('RotoDiffCWFit',1.1*r2,[],[],[],[], t, z(:,:,jw));
        err(jw)  = RotoDiffCWFit(p(jw), t, z(:,:,jw));
    end
    z0 = RotoDiffCWFit(r2, t);

    close; 
    semilogx(t,z0/z0(1,1))
    grid
    xlabel('time (ns)');
    ylabel('correlation');
    legend({'\itg\rm^{|| ||}_{|| ||}','\itg\rm_{|| ||}^{|| \perp}','\itg\rm_{|| ||}^{\perp \perp}','\itg\rm_{|| \perp}^{|| \perp}'})
    if prn
        print -depsc -r300 RotoDiffInf
    end    

    semilogx(t,[z0(:,1) squeeze(z(:,1,:))]./(ones(length(t),1)*max(abs([z0(:,1) squeeze(z(:,1,:))]))))
    colorize
    xlabel('time (ns)'); ylabel('correlation'); title('\itg\rm_{|| || \times || ||}')
    if prn
        print -depsc -r300 RotoDiffWpppp
    end
    
    semilogx(t,[z0(:,2) squeeze(z(:,2,:))]./(ones(length(t),1)*max(abs([z0(:,2) squeeze(z(:,2,:))]))))
    colorize
    xlabel('time (ns)'); ylabel('correlation'); title('\itg\rm_{|| || \times || \perp}')
    if prn
        print -depsc -r300 RotoDiffWpppo
    end
    
    semilogx(t,[z0(:,3) squeeze(z(:,3,:))]./(ones(length(t),1)*max(abs([z0(:,3) squeeze(z(:,3,:))]))))
    colorize
    xlabel('time (ns)'); ylabel('correlation'); title('\itg\rm_{|| \perp \times || ||}')
    if prn
        print -depsc -r300 RotoDiffWppop
    end
    
    semilogx(t,[z0(:,4) squeeze(z(:,4,:))]./(ones(length(t),1)*max(abs([z0(:,4) squeeze(z(:,4,:))]))))
    colorize
    xlabel('time (ns)'); ylabel('correlation'); title('\itg\rm_{|| || \times \perp \perp}')
    if prn
        print -depsc -r300 RotoDiffWppoo
    end
    
    c = polyfit(wv*1e3,100*(p-r2)/r2,2);
    plot(1e3*wv,polyval(c,wv*1e3))
    ax = axis;
    axis([1e3*min(w0) 1e3*max(w0) ax(3) ax(4)])
    grid
    xlabel('beam waist radius (nm)');
    ylabel('rel. error (%)')
    if prn
        print -depsc -r300 RotoDiffWerr
    end
end

if 0 % model curve for depolarization/focusing
    load ConfoAnisoCoefficients w0
    
    z = RotoDiffCWFit(r2, t, [], [], [], [], w0(end));
    z0 = RotoDiffCWFit(r2, t);

    close; 
    semilogx(t,z/z(1,1))
    hold on
    semilogx(t,z0/z0(1,1),'--')
    hold off
    grid
    xlabel('time (ns)');
    ylabel('correlation');
    legend({'\itg\rm^{|| ||}_{|| ||}','\itg\rm_{|| ||}^{|| \perp}','\itg\rm_{|| ||}^{\perp \perp}','\itg\rm_{|| \perp}^{|| \perp}'})
    if prn
        print -depsc -r300 RotoDiffPolarization
    end    

end

if 0 % influence of ex/em angle
    clear p
    alphav = 0:0.25:25;
    z = zeros(length(t),4,length(alphav));
    for jw=1:length(alphav)
        z(:,:,jw) = RotoDiffCWFit(r2, t, [], [], [], alphav(jw), inf);
        p(jw) = Simplex('RotoDiffCWFit',1.1*r2,[],[],[],[], t, z(:,:,jw), [], [], 0, inf);
        p(jw) = Simplex('RotoDiffCWFit',p(jw),[],[],[],[], t, z(:,:,jw), [], [], 0, inf);
        err(jw)  = RotoDiffCWFit(p(jw), t, z(:,:,jw), [], [], 0, inf);
    end
    
    close 
    
    semilogx(t,squeeze(z(:,1,:))./(ones(length(t),1)*max(abs(squeeze(z(:,1,:))))))
    colorize
    xlabel('time (ns)'); ylabel('correlation'); title('\itg\rm_{|| || \times || ||}')
    if prn
        print -depsc -r300 RotoDiffApppp
    end
    
    semilogx(t,squeeze(z(:,2,:))./(ones(length(t),1)*max(abs(squeeze(z(:,2,:))))))
    colorize
    xlabel('time (ns)'); ylabel('correlation'); title('\itg\rm_{|| || \times || \perp}')
    if prn
        print -depsc -r300 RotoDiffApppo
    end
    
    semilogx(t,squeeze(z(:,3,:))./(ones(length(t),1)*max(abs(squeeze(z(:,3,:))))))
    colorize
    xlabel('time (ns)'); ylabel('correlation'); title('\itg\rm_{|| \perp \times || ||}')
    if prn
        print -depsc -r300 RotoDiffAppop
    end
    
    semilogx(t,squeeze(z(:,4,:))./(ones(length(t),1)*max(abs(squeeze(z(:,4,:))))))
    colorize
    xlabel('time (ns)'); ylabel('correlation'); title('\itg\rm_{|| || \times \perp \perp}')
    if prn
        print -depsc -r300 RotoDiffAppoo
    end
    
    c = polyfit(alphav,100*(p-r2)/r2,5);
    plot(alphav,polyval(c,alphav))
    ax = axis;
    axis([min(alphav) max(alphav) -0.25 0.1])
    grid
    xlabel('angle between ex. and em. dipole (°)');
    ylabel('rel. error (%)')
    if prn
        print -depsc -r300 RotoDiffAerr
    end
end

if 0 % model curve for ex/em angle
    alpha = 25;
    z = RotoDiffCWFit(r2, t, [], [], [], alpha, inf);
    z0 = RotoDiffCWFit(r2, t, [], [], [], [], inf);

    close; 
    semilogx(t,z/z(1,1))
    hold on
    semilogx(t,z0/z0(1,1),'--')
    hold off
    grid
    xlabel('time (ns)');
    ylabel('correlation');
    legend({'\itg\rm^{|| ||}_{|| ||}','\itg\rm_{|| ||}^{|| \perp}','\itg\rm_{|| ||}^{\perp \perp}','\itg\rm_{|| \perp}^{|| \perp}'})
    if prn
        print -depsc -r300 RotoDiffAngle
    end    
end

if 0 % fitting anisotropic rotor isotropically
    rv = 20:80;
    z = zeros(length(t),4,length(rv));
    for jw=1:length(rv)
        z(:,:,jw) = RotoDiffCWFit([rv(jw) r2], t, [], [], [], [], inf);
        p1(jw) = Simplex('RotoDiffCWFit',1.1*r2,[],[],[],[], t, z(:,1:3,jw), [], [], 0, inf);
        p1(jw) = Simplex('RotoDiffCWFit',p1(jw),[],[],[],[], t, z(:,1:3,jw), [], [], 0, inf);
        [~,~,tmp]  = RotoDiffCWFit(p1(jw), t, z(:,1:3,jw), [], [], 0, inf);
        err1(jw) = max(max(abs(z(:,1:3,jw)-tmp)))/max(max(z(:,1:3,jw)));
        p2(jw) = Simplex('RotoDiffCWFit',1.1*r2,[],[],[],[], t, z(:,:,jw), [], [], 0, inf);
        p2(jw) = Simplex('RotoDiffCWFit',p2(jw),[],[],[],[], t, z(:,:,jw), [], [], 0, inf);
        [~,~,tmp]  = RotoDiffCWFit(p2(jw), t, z(:,:,jw), [], [], 0, inf);
        err2(jw) = max(max(abs(z(:,:,jw)-tmp)))/max(max(z(:,:,jw)));
    end
    
    close 
    
    semilogx(t,squeeze(z(:,1,:))./(ones(length(t),1)*max(abs(squeeze(z(:,1,:))))))
    colorize
    xlabel('time (ns)'); ylabel('correlation'); title('\itg\rm_{|| || \times || ||}')
    if prn
        print -depsc -r300 RotoDiffZpppp
    end
    
    semilogx(t,squeeze(z(:,2,:))./(ones(length(t),1)*max(abs(squeeze(z(:,2,:))))))
    colorize
    xlabel('time (ns)'); ylabel('correlation'); title('\itg\rm_{|| || \times || \perp}')
    if prn
        print -depsc -r300 RotoDiffZpppo
    end
    
    semilogx(t,squeeze(z(:,3,:))./(ones(length(t),1)*max(abs(squeeze(z(:,3,:))))))
    colorize
    xlabel('time (ns)'); ylabel('correlation'); title('\itg\rm_{|| \perp \times || ||}')
    if prn
        print -depsc -r300 RotoDiffZppop
    end
    
    semilogx(t,squeeze(z(:,4,:))./(ones(length(t),1)*max(abs(squeeze(z(:,4,:))))))
    colorize
    xlabel('time (ns)'); ylabel('correlation'); title('\itg\rm_{|| || \times \perp \perp}')
    if prn
        print -depsc -r300 RotoDiffZppoo
    end
    
    plot(rv,100*err1,rv,100*err2)
    grid
    legend({'CW','PIE'},2)
    xlabel('symmetry axis (Angstrom)');
    ylabel('max. rel. residue')
    if prn
        print -depsc -r300 RotoDiffAniso
    end
end

if 0 % model curves for anisotropic fitting
    z = RotoDiffCWFit([r2/2 r2], t, [], [], [], [], inf);
    p = Simplex('RotoDiffCWFit',1.1*r2,[],[],[],[], t, z, [], [], 0, inf);
    Simplex('RotoDiffCWFit',p,[],[],[],[], t, z, [], [], 0, inf);
    ax = axis;
    print -depsc -r300 RotoDiffAniso20vs40b
    p = Simplex('RotoDiffCWFit',1.1*r2,[],[],[],[], t, z(:,1:3), [], [], 0, inf);
    Simplex('RotoDiffCWFit',p,[],[],[],[], t, z(:,1:3), [], [], 0, inf);
    axis(ax);
    print -depsc -r300 RotoDiffAniso20vs40a
    
    z = RotoDiffCWFit([r2*2 r2], t, [], [], [], [], inf);
    p = Simplex('RotoDiffCWFit',1.1*r2,[],[],[],[], t, z, [], [], 0, inf);
    Simplex('RotoDiffCWFit',p,[],[],[],[], t, z, [], [], 0, inf);
    ax = axis;
    print -depsc -r300 RotoDiffAniso80vs40b
    p = Simplex('RotoDiffCWFit',1.1*r2,[],[],[],[], t, z(:,1:3), [], [], 0, inf);
    Simplex('RotoDiffCWFit',p,[],[],[],[], t, z(:,1:3), [], [], 0, inf);
    axis(ax);
    print -depsc -r300 RotoDiffAniso80vs40a
end

if 1 % experimental data: Aldolase
    clear all
    close all
    
    dirname = 'm:\MTusers\Christoph\Rotation\2011-04-05 Ald Cy5 ulong\';

    nmax = 100;
    y = zeros(241,4,nmax);
    for k=1:nmax
        for kk=1:5
            ind = round(4.5+11*rand);
            
            file = ['Point_' num2str(round(4.5+11*rand), '%d')];
            load([dirname file '.mat']);
            [y2, t] = FCSCrossReadRot(res);
            y(:,:,k) = y(:,:,k) + sum(y2,3);
        end
        k
    end

    p1inf = zeros(2,nmax);
    p2inf = zeros(2,nmax);
    p3inf = zeros(3,nmax);
    p3ww = zeros(3,nmax);
    p3full = zeros(3,nmax);
    err1inf = zeros(1,nmax);
    err2inf = err1inf;
    err3inf = err1inf;
    err3ww = err1inf;
    err3full = err1inf;

    tt = 22; % temperature
    
    for k=1:nmax
        ind = t>1e3;
        for j=1:4
            px = Simplex('ExpFun',1e3,[],[],[],[],t(ind),y(ind,j,k),1);
            [err, c] = ExpFun(px, t(ind), y(ind,j,k), 1);
            z(:,j) = y(:,j,k) - ExpFun(px, c, t);
        end
        
        ind = t<1e3;

        p1inf(:,k) = Simplex('RotoDiffCWFit',[40 200],[],[],[],[], t(ind), z(ind,1:3), [1 1], tt, 0, inf);
        for kk=1:5
            p1inf(:,k) = Simplex('RotoDiffCWFit',p1inf(:,k),[],[],[],[], t(ind), z(ind,1:3), [1 1], tt, 0, inf);
        end
        [~,~,tmp]  = RotoDiffCWFit(p1inf(:,k), t(ind), z(ind,1:3), [1 1], tt, 0, inf);
        err1inf(k) = max(max(abs(z(ind,1:3)-tmp)))/max(max(z(ind,1:3)));

        p2inf(:,k) = Simplex('RotoDiffCWFit',[40 200],[],[],[],[], t(ind), z(ind,:), [1 1], tt, 0, inf);
        for kk=1:5
            p2inf(:,k) = Simplex('RotoDiffCWFit',p2inf(:,k),[],[],[],[], t(ind), z(ind,:), [1 1], tt, 0, inf);
        end
        [~,~,tmp]  = RotoDiffCWFit(p2inf(:,k), t(ind), z(ind,:), [1 1], tt, 0, inf);
        err2inf(k) = max(max(abs(z(ind,:)-tmp)))/max(max(z(ind,:)));

        p3inf(:,k) = Simplex('RotoDiffCWFit',[40 50 200],[],[],[],[], t(ind), z(ind,1:4), [1 1], tt, 0, inf);
        for kk=1:5
            p3inf(:,k) = Simplex('RotoDiffCWFit',p3inf(:,k),[],[],[],[], t(ind), z(ind,1:4), [1 1], tt, 0, inf);
        end
        [~,~,tmp]  = RotoDiffCWFit(p3inf(:,k), t(ind), z(ind,:), [1 1], tt, 0, inf);
        err3inf(k) = max(max(abs(z(ind,:)-tmp)))/max(max(z(ind,:)));

        p3ww(:,k) = Simplex('RotoDiffCWFit',[40 50 200],[],[],[],[], t(ind), z(ind,1:4), [1 1], tt, 0, 350);
        for kk=1:5
            p3ww(:,k) = Simplex('RotoDiffCWFit',p3ww(:,k),[],[],[],[], t(ind), z(ind,1:4), [1 1], tt, 0, 350);
        end
        [~,~,tmp]  = RotoDiffCWFit(p3ww(:,k), t(ind), z(ind,:), [1 1], tt, 0, inf);
        err3ww(k) = max(max(abs(z(ind,:)-tmp)))/max(max(z(ind,:)));

        p3full(:,k) = Simplex('RotoDiffCWFit',[40 50 200],[],[],[],[], t(ind), z(ind,1:4), [1 1], tt, 10, 350);
        for kk=1:5
            p3full(:,k) = Simplex('RotoDiffCWFit',p3full(:,k),[],[],[],[], t(ind), z(ind,1:4), [1 1], tt, 10, 350);
        end
        [~,~,tmp]  = RotoDiffCWFit(p3full(:,k), t(ind), z(ind,:), [1 1], tt, 0, inf);
        err3full(k) = max(max(abs(z(ind,:)-tmp)))/max(max(z(ind,:)));
        
        disp(k)
    end
    
    save AldolaseResults2011_05_10tt22 p1inf p2inf p3inf p3ww p3full err1inf err2inf err3inf err3ww err3full y t
end

if 1 % Aldolase data evaluation

    load AldolaseResults2011_05_10tt22 
    
    mm = zeros(3,5);
    ss = mm;
    
    [h,bin] = mHist(p1inf(1,:));
    tmp = Simplex('Gauss',[mean(p1inf(1,:)) std(p1inf(1,:))],[],[],[],[],bin,h,[],[],1);
    mm(1,1) = tmp(1);
    ss(1,1) = tmp(2);
    [h,bin] = mHist(p1inf(2,:));
    tmp = Simplex('Gauss',[mean(p1inf(2,:)) std(p1inf(2,:))],[],[],[],[],bin,h,[],[],1);
    mm(3,1) = tmp(1);
    ss(3,1) = tmp(2);
    mm(2,1) = 0;
    ss(2,1) = 0;        

    [h,bin] = mHist(p2inf(1,:));
    tmp = Simplex('Gauss',[mean(p2inf(1,:)) std(p2inf(1,:))],[],[],[],[],bin,h,[],[],1);
    mm(1,2) = tmp(1);
    ss(1,2) = tmp(2);
    [h,bin] = mHist(p2inf(2,:));
    tmp = Simplex('Gauss',[mean(p2inf(2,:)) std(p2inf(2,:))],[],[],[],[],bin,h,[],[],1);
    mm(3,2) = tmp(1);
    ss(3,2) = tmp(2);
    mm(2,2) = 0;
    ss(2,2) = 0;        
   
    [h,bin] = mHist(p3inf(1,:));
    tmp = Simplex('Gauss',[mean(p3inf(1,:)) std(p3inf(1,:))],[],[],[],[],bin,h,[],[],1);
    mm(1,3) = tmp(1);
    ss(1,3) = tmp(2);
    [h,bin] = mHist(p3inf(2,:));
    tmp = Simplex('Gauss',[mean(p3inf(2,:)) std(p3inf(2,:))],[],[],[],[],bin,h,[],[],1);
    mm(2,3) = tmp(1);
    ss(2,3) = tmp(2);
    [h,bin] = mHist(p3inf(3,:));
    tmp = Simplex('Gauss',[mean(p3inf(3,:)) std(p3inf(3,:))],[],[],[],[],bin,h,[],[],1);
    mm(3,3) = tmp(1);
    ss(3,3) = tmp(2);        
    
    [h,bin] = mHist(p3ww(1,:));
    tmp = Simplex('Gauss',[mean(p3ww(1,:)) std(p3ww(1,:))],[],[],[],[],bin,h,[],[],1);
    mm(1,4) = tmp(1);
    ss(1,4) = tmp(2);
    [h,bin] = mHist(p3ww(2,:));
    tmp = Simplex('Gauss',[mean(p3ww(2,:)) std(p3ww(2,:))],[],[],[],[],bin,h,[],[],1);
    mm(2,4) = tmp(1);
    ss(2,4) = tmp(2);
    [h,bin] = mHist(p3ww(3,:));
    tmp = Simplex('Gauss',[mean(p3ww(3,:)) std(p3ww(3,:))],[],[],[],[],bin,h,[],[],1);
    mm(3,4) = tmp(1);
    ss(3,4) = tmp(2);        
    
    [h,bin] = mHist(p3full(1,:));
    tmp = Simplex('Gauss',[mean(p3full(1,:)) std(p3full(1,:))],[],[],[],[],bin,h,[],[],1);
    mm(1,5) = tmp(1);
    ss(1,5) = tmp(2);
    [h,bin] = mHist(p3full(2,:));
    tmp = Simplex('Gauss',[mean(p3full(2,:)) std(p3full(2,:))],[],[],[],[],bin,h,[],[],1);
    mm(2,5) = tmp(1);
    ss(2,5) = tmp(2);
    [h,bin] = mHist(p3full(3,:));
    tmp = Simplex('Gauss',[mean(p3full(3,:)) std(p3full(3,:))],[],[],[],[],bin,h,[],[],1);
    mm(3,5) = tmp(1);
    ss(3,5) = tmp(2);    
    
end