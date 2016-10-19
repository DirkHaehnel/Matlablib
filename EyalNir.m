% program Eyal Nir

clear all
close all
zpin = 0;
atf = [];
kappa = 1;
lt = 0;
sat = 0;

if 0 % 12.6.
    NA = 1.2;
    ng = 1.52;
    nw = 1.333;
    mag = 60;
    over = [NA*180/mag*1e3 0 3e3];
    focpos = 20;
    av = 100/2;
    lamexv = [0.532 0.647];
    lamemv = [0.550 0.650];
    dv = 0:5;
    
    for jlam=1:length(lamexv)
        lamex = lamexv(jlam);
        lamem = lamemv(jlam);        
        for jd = 1:length(dv)
            atf = [ng dv(jd)];
            zfield = [15 25]; 
            rhofield = [0 5];    
            resolution = 30;
            [volx, voly, rho, z, h0, h1, h2] = ...
                MDF(rhofield, zfield, NA, nw, nw, nw, [], 0, [], lamex, over, focpos, lamem, mag, av, zpin, atf, kappa, lt, sat, resolution);
            [veffc(jd,jlam), veff(jd,jlam)] = DetectionVolume(rho,z,volx,voly);
            feld = [volx(:,:,1)+voly(:,:,1)];
            mpc(rho,z,feld,'lr');
            xlabel('\rho [\mum]');
            ylabel('\itz\rm [\mum]');
            title(['cov. slide ' mint2str(dv(jd),1) ' \mum: V = ' mnum2str(veff(jd,jlam),1,3) ' fl']);
            eval(['print -dpng EyalNirA' mint2str(1e3*lamex,3) 'ex' mint2str(1e3*lamem,3) 'em' mint2str(dv(jd),1) 'd'])
        end
    end
    dd=0:0.01:5;
    plot(dd,interp1(dv,veff,dd,'cubic'))
    xlabel('cover-slide thickness deviation [\mum]')
    ylabel('detection volume [\mum^3]')
    print -dpng EyalNirA 
end

if 0 % 12.6.
    NA = 1.2;
    ng = 1.52;
    nw = 1.333;
    mag = 60;
    focpos = 20;
    av = 100/2;
    lamexv = [0.532 0.647];
    lamemv = [0.550 0.650];
    dv = 0:5;
    laser = [2e3 4e3]; 
    
    for jlam=1:length(lamexv)
        lamex = lamexv(jlam);
        lamem = lamemv(jlam);        
        over = [NA*180/mag*1e3 0 laser(jlam)];
        for jd = 1:length(dv)
            atf = [ng dv(jd)];
            zfield = [15 25]; 
            rhofield = [0 5];    
            resolution = 30;
            [volx, voly, rho, z, h0, h1, h2] = ...
                MDF(rhofield, zfield, NA, nw, nw, nw, [], 0, [], lamex, over, focpos, lamem, mag, av, zpin, atf, kappa, lt, sat, resolution);
            [veffc(jd,jlam), veff(jd,jlam)] = DetectionVolume(rho,z,volx,voly);
            feld = [volx(:,:,1)+voly(:,:,1)];
            mpc(rho,z,feld,'lr');
            xlabel('\rho [\mum]');
            ylabel('\itz\rm [\mum]');
            title(['cov. slide ' mint2str(dv(jd),1) ' \mum: V = ' mnum2str(veff(jd,jlam),1,3) ' fl']);
            eval(['print -dpng EyalNirB' mint2str(1e3*lamex,3) 'ex' mint2str(1e3*lamem,3) 'em' mint2str(dv(jd),1) 'd'])
        end
    end
    dd=0:0.01:5;
    plot(dd,interp1(dv,veff,dd,'cubic'))
    xlabel('cover-slide thickness deviation [\mum]')
    ylabel('detection volume [\mum^3]')
    print -dpng EyalNirB 
end

if 0 % 12.6.
    NA = 1.2;
    ng = 1.52;
    nw = 1.333;
    mag = 60;
    focpos = 20;
    avv = [100 80]/2;
    lamexv = [0.532 0.647];
    lamemv = [0.550 0.650];
    dv = 0:5;
    over = [NA*180/mag*1e3 0 3e3];
    
    for jlam=1:length(lamexv)
        lamex = lamexv(jlam);
        lamem = lamemv(jlam);
        av = avv(jlam);
        for jd = 1:length(dv)
            atf = [ng dv(jd)];
            zfield = [15 25]; 
            rhofield = [0 5];    
            resolution = 30;
            [volx, voly, rho, z, h0, h1, h2] = ...
                MDF(rhofield, zfield, NA, nw, nw, nw, [], 0, [], lamex, over, focpos, lamem, mag, av, zpin, atf, kappa, lt, sat, resolution);
            [veffc(jd,jlam), veff(jd,jlam)] = DetectionVolume(rho,z,volx,voly);
            feld = [volx(:,:,1)+voly(:,:,1)];
            mpc(rho,z,feld,'lr');
            xlabel('\rho [\mum]');
            ylabel('\itz\rm [\mum]');
            title(['cov. slide ' mint2str(dv(jd),1) ' \mum: V = ' mnum2str(veff(jd,jlam),1,3) ' fl']);
            eval(['print -dpng EyalNirC' mint2str(1e3*lamex,3) 'ex' mint2str(1e3*lamem,3) 'em' mint2str(dv(jd),1) 'd'])
        end
    end
    dd=0:0.01:5;
    plot(dd,interp1(dv,veff,dd,'cubic'))
    xlabel('cover-slide thickness deviation [\mum]')
    ylabel('detection volume [\mum^3]')
    print -dpng EyalNirC 
end

cd d:/joerg/doc/eyalnir/daten
if 0 % 19.6. / FRET Cy3b and Atto647; pinhole size and cover slide
    NA = 1.2;
    ng = 1.52;
    nw = 1.333;
    mag = 60;
    focpos = 20;
    avv = (50:10:150)/2;
    lamex = 0.532;
    lamgreen = 0.570;
    lamred = 0.670;
    dv = 0:5;
    zpinv = [-5e3:1e3:5e3];
    laser = 3e3; 
    over = [NA*180/mag*1e3 0 laser];
    zfield = [15 25]; 
    rhofield = [0 5];    
    resolution = 30;

    jzpin0 = 1;
    jd0 = 1;
    jav0 = 1;
    
    for jzpin = max([1 jzpin0]):length(zpinv)
        jzpin0 = -1;
        zpin = zpinv(jzpin);
        for jd = max([1 jd0]):length(dv)
            jd0 = -1;
            d = dv(jd);
            atf = [ng d];
            for jav=max([1 jav0]):length(avv)
                jav0 = -1;
                av = avv(jav);            
                lamem = lamgreen;  
                [volxg, volyg, rho, z] = ...
                    MDF(rhofield, zfield, NA, nw, nw, nw, [], 0, [], lamex, over, focpos, lamem, mag, av, zpin, atf, kappa, lt, sat, resolution);
                if zpin<0
                    eval(['save EyalNirGreen' mint2str(av*2,3) 'av' mint2str(d,1) 'd-' mint2str(zpin/1e3,1) 'zpin rho z volxg volyg']);
                else
                    eval(['save EyalNirGreen' mint2str(av*2,3) 'av' mint2str(d,1) 'd+' mint2str(zpin/1e3,1) 'zpin rho z volxg volyg']);
                end
                [tmp, vgreen(jav,jd,jzpin)] = DetectionVolume(rho,z,volxg,volyg);
                lamem = lamred;
                [volxr, volyr] = ...
                    MDF(rhofield, zfield, NA, nw, nw, nw, [], 0, [], lamex, over, focpos, lamem, mag, av, zpin, atf, kappa, lt, sat, resolution);
                if zpin<0
                    eval(['save EyalNirRed' mint2str(av*2,3) 'av' mint2str(d,1) 'd-' mint2str(zpin/1e3,1) 'zpin rho z volxr volyr']);
                else
                    eval(['save EyalNirRed' mint2str(av*2,3) 'av' mint2str(d,1) 'd+' mint2str(zpin/1e3,1) 'zpin rho z volxr volyr']);
                end
            end
        end
    end
end


if 1 % 27.6. / FRET Cy3b and Atto647; pinhole size and cover slide
    NA = 1.2;
    ng = 1.52;
    nw = 1.333;
    mag = 60;
    focpos = 20;
    avv = (50:10:150)/2;
    lamex = 0.532;
    lamgreen = 0.570;
    lamred = 0.670;
    dv = 0:5;
    zpinv = [-5e3:1e3:5e3];
    laser = 3e3; 
    over = [NA*180/mag*1e3 0 laser];
    zfield = [15 25]; 
    rhofield = [0 5];    
    resolution = 30;

    load EyalNirStart
    jzpin0 = jzpin;
    jd0 = jd;
    jav0 = jav;
    kzpin0 = kzpin;
    kd0 = kd;
    kav0 = kav;
    
    for jzpin = max([1 jzpin0]):length(zpinv)
        jzpin0 = -1;
        zping = zpinv(jzpin);
        for jd = max([1 jd0]):length(dv)
            jd0 = -1;
            dg = dv(jd);
            atfg = [ng dg];
            for jav=max([1 jav0]):length(avv)
                jav0 = -1;
                avg = avv(jav);
                if zping<0
                    eval(['load EyalNirGreen' mint2str(avg*2,3) 'av' mint2str(dg,1) 'd-' mint2str(zping/1e3,1) 'zpin rho z volxg volyg']);
                else
                    eval(['load EyalNirGreen' mint2str(avg*2,3) 'av' mint2str(dg,1) 'd+' mint2str(zping/1e3,1) 'zpin rho z volxg volyg']);
                end
                [tmp, vgreen(jav,jd,jzpin)] = DetectionVolume(rho,z,volxg,volyg);
                if zping<0
                    eval(['load EyalNirRed' mint2str(avg*2,3) 'av' mint2str(dg,1) 'd-' mint2str(zping/1e3,1) 'zpin rho z volxr volyr']);
                else
                    eval(['load EyalNirRed' mint2str(avg*2,3) 'av' mint2str(dg,1) 'd+' mint2str(zping/1e3,1) 'zpin rho z volxr volyr']);
                end
                [tmp, vred(jav,jd,jzpin)] = DetectionVolume(rho,z,volxr,volyr);
                for kzpin = max([1 kzpin0]):length(zpinv)
                    kzpin0 = -1;
                    zpinr = zpinv(kzpin);
                    for kd = max([1 kd0]):length(dv)
                        kd0 = -1;
                        dr = dv(kd);
                        atfr = [ng dr];
                        for kav=max([1 kav0]):length(avv)
                            kav0 = -1;
                            avr = avv(kav);
                            if zpinr<0
                                eval(['load EyalNirRed' mint2str(avr*2,3) 'av' mint2str(dr,1) 'd-' mint2str(zpinr/1e3,1) 'zpin rho z volxr volyr']);
                            else
                                eval(['load EyalNirRed' mint2str(avr*2,3) 'av' mint2str(dr,1) 'd+' mint2str(zpinr/1e3,1) 'zpin rho z volxr volyr']);
                            end
                            vcross(jav,jd,jzpin,kav,kd,kzpin) = DetectionVolume(rho,z,volxg+volyg,volxr+volyr);
                        end
                    end
                end
                save EyalNir vgreen vred vcross
            end
        end
    end
end

