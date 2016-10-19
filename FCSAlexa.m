cd D:\Joerg\Doc\Fcs\OptSat

if 1
    close all
    clear all
    NA = 1.2;
    wd = 3e3;
    tubelens = 180e3;
    n0 = 1.333;
    n = 1.333;
    n1 = 1.333;
    d0 = [];
    d = [];
    d1 = [];
    lamex = 0.635;
    focpos = 10;
    lamem = 0.67;
    mag = tubelens/wd; %60;
    av = 100/2; %pinhole radius!!!
    zpin = 0; 
    atf = [1.52 0];
    kappa = 0; 
    lt = [];
    sat = 0:0.1:15;

    ww = 4.9e3;
    for j = 1:5
        zeta = (j-1)/20*2*pi*ww^2/NA^2/wd^2;
        w0 = ww/sqrt(1+zeta^2);
        z0 = pi*w0^2*zeta/lamex;
        over = [wd*NA z0 -z0 w0 w0];
        
        [modres(:,:,j), autotime, veff(1,:,j), exc(:,:,j)] = NewFCS(NA, n0, n, n1, d0, d, d1, lamex, lamem, mag, focpos, av, over, zpin, kappa, atf, lt, sat);
        save AlexaModel49_0
    end

    clear modres veff exc
    ww = 5.3e3;
    for j = 1:5
        zeta = (j-1)/20*2*pi*ww^2/NA^2/wd^2;
        w0 = ww/sqrt(1+zeta^2);
        z0 = pi*w0^2*zeta/lamex;
        over = [wd*NA z0 -z0 w0 w0];
        
        [modres(:,:,j), autotime, veff(1,:,j), exc(:,:,j)] = NewFCS(NA, n0, n, n1, d0, d, d1, lamex, lamem, mag, focpos, av, over, zpin, kappa, atf, lt, sat);
        save AlexaModel53_0
    end
    
    
% 
%     ww = 5.3e3;
%     for j = 1:5
%         atf(2) = 2;
%         zeta = (j-1)/20*2*pi*ww^2/NA^2/wd^2;
%         w0 = ww/sqrt(1+zeta^2);
%         z0 = pi*w0^2*zeta/lamex;
%         over = [wd*NA z0 -z0 w0 w0];
%         
%         [modres(:,:,j), autotime, veff(1,:,j), exc(:,:,j)] = NewFCS(NA, n0, n, n1, d0, d, d1, lamex, lamem, mag, focpos, av, over, zpin, kappa, atf, lt, sat);
%     end
%     save pulsemodel53_2
% 
%     ww = 5.3e3;
%     for j = 1:5
%         atf(2) = 4;
%         zeta = (j-1)/20*2*pi*ww^2/NA^2/wd^2;
%         w0 = ww/sqrt(1+zeta^2);
%         z0 = pi*w0^2*zeta/lamex;
%         over = [wd*NA z0 -z0 w0 w0];
%         
%         [modres(:,:,j), autotime, veff(1,:,j), exc(:,:,j)] = NewFCS(NA, n0, n, n1, d0, d, d1, lamex, lamem, mag, focpos, av, over, zpin, kappa, atf, lt, sat);
%     end
%     save pulsemodel53_4
% 
%     ww = 5.3e3;
%     for j = 1:5
%         atf(2) = 6;
%         zeta = (j-1)/20*2*pi*ww^2/NA^2/wd^2;
%         w0 = ww/sqrt(1+zeta^2);
%         z0 = pi*w0^2*zeta/lamex;
%         over = [wd*NA z0 -z0 w0 w0];
%         
%         [modres(:,:,j), autotime, veff(1,:,j), exc(:,:,j)] = NewFCS(NA, n0, n, n1, d0, d, d1, lamex, lamem, mag, focpos, av, over, zpin, kappa, atf, lt, sat);
%     end
%     save pulsemodel53_6

%     ww = 5.3e3;
%     for j = 1:5
%         atf(2) = 8;
%         zeta = (j-1)/20*2*pi*ww^2/NA^2/wd^2;
%         w0 = ww/sqrt(1+zeta^2);
%         z0 = pi*w0^2*zeta/lamex;
%         over = [wd*NA z0 -z0 w0 w0];
%         
%         [modres(:,:,j), autotime, veff(1,:,j), exc(:,:,j)] = NewFCS(NA, n0, n, n1, d0, d, d1, lamex, lamem, mag, focpos, av, over, zpin, kappa, atf, lt, sat);
%     end
%     save pulsemodel53_8
%     
%     ww = 5.3e3;
%     for j = 1:5
%         atf(2) = 10;
%         zeta = (j-1)/20*2*pi*ww^2/NA^2/wd^2;
%         w0 = ww/sqrt(1+zeta^2);
%         z0 = pi*w0^2*zeta/lamex;
%         over = [wd*NA z0 -z0 w0 w0];
%         
%         [modres(:,:,j), autotime, veff(1,:,j), exc(:,:,j)] = NewFCS(NA, n0, n, n1, d0, d, d1, lamex, lamem, mag, focpos, av, over, zpin, kappa, atf, lt, sat);
%     end
%     save pulsemodel53_10
    
end

if 0
    close all
    clear all
    NA = 1.2;
    wd = 3e3;
    tubelens = 180e3;
    n0 = 1.333;
    n = 1.333;
    n1 = 1.333;
    d0 = [];
    d = [];
    d1 = [];
    lamex = 0.635;
    focpos = 10;
    lamem = 0.67;
    mag = tubelens/wd; %60;
    av = 100/2; %pinhole radius!!!
    zpin = 0; 
    atf = [1.51 0];
    kappa = 0; 
    lt = [];
    
    sat = 0:0.1:6;
    w0 = [500.9517 841.7505];
    z0 = 1e7*[1.1879 2.0499];
    over = [wd*NA z0 w0];
    [modres, autotime, veff, exc, rho, z, volx, voly] = NewFCS(NA, n0, n, n1, d0, d, d1, lamex, lamem, mag, focpos, av, over, zpin, kappa, atf, lt, sat);
    save pulsemodel
    
    sat = 0:0.1:6;
    w0 = [3.1514    1.0195]*1e3;
    z0 = [4.4967   -1.9257]*1e7;
    over = [wd*NA z0([2 1]) w0([2 1])];
    [modres, autotime, veff, exc, rho, z, volx, voly] = NewFCS(NA, n0, n, n1, d0, d, d1, lamex, lamem, mag, focpos, av, over, zpin, kappa, atf, lt, sat);
    save cwmodel2
end


% plot(0:0.1:1,exc(3:4,:)./(ones(2,1)*sum(exc(3:4,:))))
% plot(0:0.1:1,(exc(4,:)+exc(5,:))./exc(3,:),0:0.1:1,(exc(6,:)-(exc(4,:)+exc(5,:)).^2)./exc(3,:))

if 0
    clear all
    load F:\Daten\Gregor\040316\alexapower.mat 
    load cwmodeltotal
    
    tt = 10:length(t);
    for jm = 1:size(modres,3)
        for j=1:size(x,2)
            for js=1:length(sat)
                p(:,js,j,jm) = simplex('alexafit',[5e-7 1e5],zeros(1,2),[],[],[],exc(:,js),t(tt),x(tt,j),autotime,squeeze(modres(:,js,jm)));
                p(:,js,j,jm) = simplex('alexafit',p(:,js,j,jm),zeros(1,2),[],[],[],exc(:,js),t(tt),x(tt,j),autotime,squeeze(modres(:,js,jm)));
                [err(js,j,jm), cc(:,js,j,jm), y(:,js,j,jm)] = alexafit(p(:,js,j,jm),exc(:,js),t(tt),x(tt,j),autotime,squeeze(modres(:,js,jm)));
            end
            
            tst = err(1,j,jm);
            ind(j,jm) = length(sat);
            for js=1:length(sat)
                if err(js,j,jm)<tst
                    tst = err(js,j,jm);
                    ind(j,jm) = js;
                end
            end
        end
    end

    save cwalexatotalres
end
    
if 0
    clear all
    load F:\Daten\Gregor\040316\alexapower.mat 
    load cwmodel2

    tt = 10:length(t);    
    for j=1:size(x,2)
        for js=1:length(sat)
            p(:,js,j) = simplex('alexafit',[5e-7 1e5],zeros(1,2),[],[],[],exc(:,js),t(tt),x(tt,j),autotime,squeeze(modres(:,js)));
            p(:,js,j) = simplex('alexafit',p(:,js,j),zeros(1,2),[],[],[],exc(:,js),t(tt),x(tt,j),autotime,squeeze(modres(:,js)));
            [err(js,j), cc(:,js,j), y(:,js,j)] = alexafit(p(:,js,j),exc(:,js),t(tt),x(tt,j),autotime,squeeze(modres(:,js)));
        end
        
        tst = err(1,j);
        ind(j) = length(sat);
        for js=1:length(sat)
            if err(js,j)<tst
                tst = err(js,j);
                ind(j) = js;
            end
        end
    end
    
    save cwalexa2res
    
    load cwmodel1

    tt = 10:length(t);    
    for j=1:size(x,2)
        for js=1:length(sat)
            p(:,js,j) = simplex('alexafit',[5e-7 1e5],zeros(1,2),[],[],[],exc(:,js),t(tt),x(tt,j),autotime,squeeze(modres(:,js)));
            p(:,js,j) = simplex('alexafit',p(:,js,j),zeros(1,2),[],[],[],exc(:,js),t(tt),x(tt,j),autotime,squeeze(modres(:,js)));
            [err(js,j), cc(:,js,j), y(:,js,j)] = alexafit(p(:,js,j),exc(:,js),t(tt),x(tt,j),autotime,squeeze(modres(:,js)));
        end
        
        tst = err(1,j);
        ind(j) = length(sat);
        for js=1:length(sat)
            if err(js,j)<tst
                tst = err(js,j);
                ind(j) = js;
            end
        end
    end
    
    save cwalexa1res
    
end

% jm=5; for j=1:length(sat) tst(j,3)=5.5e-13./p(1,j,5,jm); tst(j,1)=interp1(sat,5.5e-13./p(1,:,3,jm),sat(j)/10,'cubic'); tst(j,2)=interp1(sat,5.5e-13./p(1,:,4,jm),3*sat(j)/10,'cubic'); end

if 0
    clear all
    load F:\Daten\Gregor\040309\alexa.mat
    load pulsemodeltotal
    
    tt = 10:length(t);
    for jm = 1:size(modres,3)
        for j=1:size(x,2)
            for js=1:length(sat)
                p(:,js,j,jm) = simplex('alexafit',[5e-7 1e5],zeros(1,2),[],[],[],exc(:,js),t(tt),x(tt,j),autotime,squeeze(modres(:,js,jm)));
                p(:,js,j,jm) = simplex('alexafit',p(:,js,j,jm),zeros(1,2),[],[],[],exc(:,js),t(tt),x(tt,j),autotime,squeeze(modres(:,js,jm)));
                [err(js,j,jm), cc(:,js,j,jm), y(:,js,j,jm)] = alexafit(p(:,js,j,jm),exc(:,js),t(tt),x(tt,j),autotime,squeeze(modres(:,js,jm)));
            end
            
            tst = err(1,j,jm);
            ind(j,jm) = length(sat);
            for js=1:length(sat)
                if err(js,j,jm)<tst
                    tst = err(js,j,jm);
                    ind(j,jm) = js;
                end
            end
        end
    end

    save pulsealexatotalres2
end

if 0
    tt = 10:length(t);
    for jm = 1:size(modres,3)
        for j=1:size(x,2)
            for js=1:length(sat)
                p(:,js,j,jm) = simplex('affinefit',5e-7,0,[],[],[],t(tt),x(tt,j),autotime,squeeze(modres(:,js,jm)));
                p(:,js,j,jm) = simplex('affinefit',p(:,js,j,jm),0,[],[],[],t(tt),x(tt,j),autotime,squeeze(modres(:,js,jm)));
                [err(js,j,jm), cc(:,js,j,jm), y(:,js,j,jm)] = affinefit(p(:,js,j,jm),t(tt),x(tt,j),autotime,squeeze(modres(:,js,jm)));
            end
            
            tst = err(1,j,jm);
            ind(j,jm) = length(sat);
            for js=1:length(sat)
                if err(js,j,jm)<tst
                    tst = err(js,j,jm);
                    ind(j,jm) = js;
                end
            end
        end
    end
end

if 0
    load F:\Daten\Gregor\040309\bla.mat
    
    for jatf = 0:2:6
        load(['pulsemodel53_' num2str(jatf)]);
        
        tt = 10:length(t);
        for jm = 1:size(modres,3)
            for j=1:size(x,2)
                for js=1:length(sat)
                    p(:,js,j,jm) = simplex('affineexpfit',[5e-7 1e5],[0 0],[],[],[],t(tt),x(tt,j),autotime,squeeze(modres(:,js,jm)));
                    p(:,js,j,jm) = simplex('affineexpfit',p(:,js,j,jm),[0 0],[],[],[],t(tt),x(tt,j),autotime,squeeze(modres(:,js,jm)));
                    [err(js,j,jm), cc(:,js,j,jm), y(:,js,j,jm)] = affineexpfit(p(:,js,j,jm),t(tt),x(tt,j),autotime,squeeze(modres(:,js,jm)));
                end
                
                tst = err(1,j,jm);
                ind(j,jm) = length(sat);
                for js=1:length(sat)
                    if err(js,j,jm)<tst
                        tst = err(js,j,jm);
                        ind(j,jm) = js;
                    end
                end
            end
        end
    
        eval(['save blares' num2str(jatf) ' p err cc y tst ind sat x t tt']);
    end
    
    load F:\Daten\Gregor\040309\lys.mat
    
    for jatf = 0:2:6
        load(['pulsemodel53_' num2str(jatf)]);
        
        tt = 10:length(t);
        for jm = 1:size(modres,3)
            for j=1:size(x,2)
                for js=1:length(sat)
                    p(:,js,j,jm) = simplex('affineexpfit',[5e-7 1e5],[0 0],[],[],[],t(tt),x(tt,j),autotime,squeeze(modres(:,js,jm)));
                    p(:,js,j,jm) = simplex('affineexpfit',p(:,js,j,jm),[0 0],[],[],[],t(tt),x(tt,j),autotime,squeeze(modres(:,js,jm)));
                    [err(js,j,jm), cc(:,js,j,jm), y(:,js,j,jm)] = affineexpfit(p(:,js,j,jm),t(tt),x(tt,j),autotime,squeeze(modres(:,js,jm)));
                end
                
                tst = err(1,j,jm);
                ind(j,jm) = length(sat);
                for js=1:length(sat)
                    if err(js,j,jm)<tst
                        tst = err(js,j,jm);
                        ind(j,jm) = js;
                    end
                end
            end
        end
    
        eval(['save lysres' num2str(jatf) ' p err cc y tst ind sat x t tt']);
    end
    
end