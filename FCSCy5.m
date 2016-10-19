close all
clear all
nn = 3;
NA = 1.2;
wd = 3e3;
tubelens = 180e3;
n0 = 1.333;
n = 1.333;
n1 = 1.333;
d0 = [];
d = 0;
d1 = [];
lamex = 0.635;
focpos = 10;
lamem = 0.67;
mag = tubelens/wd; %60;
av = 100/2; %pinhole radius!!!

kappa = 1; 
lt = [];
satv = [0 0.01*exp(log(1e3)/10*(0:10))];

cd D:\Daten\Cy5DynamicsTheory

zpin = 0; 

if 0
    for jww=1:3
        ww = 3.5e3 + 500*(jww-1);
        for jzeta = 1:5
            zeta = (jzeta-1)/20*2*pi*ww^2/NA^2/wd^2;
            w0 = ww/sqrt(1+zeta^2);
            z0 = pi*w0^2*zeta/lamex;
            over = [wd*NA z0 -z0 w0 w0];
            for jatf = 1:21
                atf = [1.52 (jatf-1)];
                for s=1:length(satv)
                    sat = satv(s);
                    outfile = ['Cy5Model' mint2str(ww/100,2) 'ww' mint2str((jzeta-1),2) 'zeta' mint2str((jatf-1),2) 'atf' mint2str(s,2) 'sat'];
                    exist([outfile '.mat'])
                    if ~(exist([outfile '.mat'])==2)
                        [volx, voly, rho, z, h0, h1, h2] = ...
                            MDF(0, [focpos-min([focpos,10]) focpos+10], NA, n0, n, n1, d0, d, d1, lamex, over, focpos, lamem, mag, av, zpin, atf, kappa, lt, sat, 10);
                        exc = volx(:,:,1,1)+voly(:,:,1,1); exc = exc/max(exc);
                        z(exc<1e-3) = [];
                        zfield = [z(1) z(end)];
                        rhofield = [0 diff(zfield)/3];
                        [volx, voly, rho, z, h0, h1, h2] = ...
                            MDF(rhofield, zfield, NA, n0, n, n1, d0, d, d1, lamex, over, focpos, lamem, mag, av, zpin, atf, kappa, lt, sat);
                        [modres, autotime] = FCS(rho, z, volx, voly, h0, h1, h2, [], false);
                        eval(['save ' outfile ' volx voly rho z h0 h1 h2 modres autotime NA wd tubelens n0 n n1 d0 d d1 lamex focpos lamem mag av zpin atf kappa lt sat'])
                    end
                end    
            end
        end
    end
end


if 1
    load D:\Joerg\Doc\Gregor\Cy5Dynamics\cy5dynamics.mat
    % t ~ autotime;
    % x ~ auto;

    for jww=1:3
        ww = 3.5e3 + 500*(jww-1);
        for jzeta = 1:5
            zeta = (jzeta-1)/20*2*pi*ww^2/NA^2/wd^2;
            w0 = ww/sqrt(1+zeta^2);
            z0 = pi*w0^2*zeta/lamex;
            over = [wd*NA z0 -z0 w0 w0];
            for jatf = 1:21
                atf = [1.52 (jatf-1)];            
                for s=1:length(satv) 
                    sat = satv(s);
                    eval(['load D:\Daten\Cy5DynamicsTheory\Cy5Model' mint2str(ww/100,2) 'ww' mint2str((jzeta-1),2) 'zeta' mint2str((jatf-1),2) ...
                            'atf' mint2str(s,2) 'sat autotime modres']);
                    for j=1:size(x,2)
                        p(:,j,s,jatf,jzeta,jww) = simplex('affineexpfit',[1e-7 1e5],[0 0],[],[],[],t,x(:,j),autotime,modres);
                        p(:,j,s,jatf,jzeta,jww) = simplex('affineexpfit',p(:,j,s,jatf,jzeta,jww),[0 0],[],[],[],t,x(:,j),autotime,modres);
                        err(j,s,jatf,jzeta,jww) = affineexpfit(p(:,j,s,jatf,jzeta,jww),t,x(:,j),autotime,modres);
                    end
                end    
            end
            save Cy5FitResults p err
        end
    end
end

if 1
    pow = [10 25 50 100 250 500]/500;
    % load D:\Daten\Cy5DynamicsTheory\Cy5Fit-4zpinResults.mat
    load D:\Daten\Cy5DynamicsTheory\Cy5FitResults.mat
    tst = 1;
    for jww=1:3
        ww = 3.5e3 + 500*(jww-1);
        for jzeta = 1:5
            zeta = (jzeta-1)/20*2*pi*ww^2/NA^2/wd^2;
            w0 = ww/sqrt(1+zeta^2);
            z0 = pi*w0^2*zeta/lamex;
            over = [wd*NA z0 -z0 w0 w0];
            for jatf = 1:21
                atf = [1.52 (jatf-1)];            
                for fac=0.5:0.1:20
                    tmp = interp1(satv,err(1,:,jatf,jzeta,jww),fac*pow(1),'cubic'); 
                    for j=2:6
                        tmp = tmp*interp1(satv,err(j,:,jatf,jzeta,jww),fac*pow(j),'cubic');
                    end
                    if tmp<tst
                        tst = tmp;
                        ind = [fac jatf jzeta jww];
                    end
                end    
            end
        end
        jww
    end
end

if 0
    pow = [10 25 50 100 250 500]/500;
    load D:\Daten\Cy5DynamicsTheory\Cy5FitResults.mat
    tstall = ones(1,6);
    for jww=1:3
        ww = 3.5e3 + 500*(jww-1);
        for jzeta = 1:5
            zeta = (jzeta-1)/20*2*pi*ww^2/NA^2/wd^2;
            w0 = ww/sqrt(1+zeta^2);
            z0 = pi*w0^2*zeta/lamex;
            over = [wd*NA z0 -z0 w0 w0];
            for jatf = 1:21
                atf = [1.52 (jatf-1)];
                for fac=0.5:0.1:20
                    for j=1:6
                        tmp = interp1(satv,err(j,:,jatf,jzeta,jww),fac*pow(j),'cubic');
                        if tmp<tstall(j)
                            tstall(j) = tmp;
                            indall(j,:) = [fac jatf jzeta jww];
                        end

                    end
                end
            end
        end
        jww
    end
end

if 0
    pow = [10 25 50 100 250 500]/500;
    load D:\Daten\Cy5DynamicsTheory\Cy5FitResults.mat
    tst5 = 1;
    for jww=1:3
        ww = 3.5e3 + 500*(jww-1);
        for jzeta = 1:5
            zeta = (jzeta-1)/20*2*pi*ww^2/NA^2/wd^2;
            w0 = ww/sqrt(1+zeta^2);
            z0 = pi*w0^2*zeta/lamex;
            over = [wd*NA z0 -z0 w0 w0];
            for jatf = 1:21
                atf = [1.52 (jatf-1)];
                for fac=0.5:0.1:20
                    tmp = interp1(satv,err(1,:,jatf,jzeta,jww),fac*pow(1),'cubic'); 
                    for j=2:5
                        tmp = tmp*interp1(satv,err(j,:,jatf,jzeta,jww),fac*pow(j),'cubic');
                    end
                    if tmp<tst5
                        tst5 = tmp;
                        ind5 = [fac jatf jzeta jww];
                    end
                end
            end
        end
        jww
    end
    tst4 = 1;
    for jww=1:3
        ww = 3.5e3 + 500*(jww-1);
        for jzeta = 1:5
            zeta = (jzeta-1)/20*2*pi*ww^2/NA^2/wd^2;
            w0 = ww/sqrt(1+zeta^2);
            z0 = pi*w0^2*zeta/lamex;
            over = [wd*NA z0 -z0 w0 w0];
            for jatf = 1:21
                atf = [1.52 (jatf-1)];
                for fac=0.5:0.1:20
                    tmp = interp1(satv,err(1,:,jatf,jzeta,jww),fac*pow(1),'cubic'); 
                    for j=2:4
                        tmp = tmp*interp1(satv,err(j,:,jatf,jzeta,jww),fac*pow(j),'cubic');
                    end
                    if tmp<tst4
                        tst4 = tmp;
                        ind4 = [fac jatf jzeta jww];
                    end
                end
            end
        end
        jww
    end
end

if 0
    tmp = interp1(satv,squeeze(p(1,1,:,ind(2),ind(3),ind(4))),fac*pow(1),'cubic');
    for j=2:6
        tmp = [tmp interp1(satv,squeeze(p(1,j,:,ind(2),ind(3),ind(4))),fac*pow(j),'cubic')];
    end
end

if 0
    pow = [10 25 50 100 250 500]/500;
    load D:\Daten\Cy5DynamicsTheory\Cy5FitResults.mat
    tst = 1;
    for jww=1:3
        ww = 3.5e3 + 500*(jww-1);
        for jzeta = 1:5
            zeta = (jzeta-1)/20*2*pi*ww^2/NA^2/wd^2;
            w0 = ww/sqrt(1+zeta^2);
            z0 = pi*w0^2*zeta/lamex;
            over = [wd*NA z0 -z0 w0 w0];
            for jatf = 1:11
                atf = [1.52 (jatf-1)];            
                for fac=0.5:0.1:10
                    tmp(1) = interp1(satv,squeeze(p(1,1,:,jatf,jzeta,jww)),fac*pow(1),'cubic'); 
                    for j=2:6
                        tmp(j) = interp1(satv,squeeze(p(1,j,:,jatf,jzeta,jww)),fac*pow(j),'cubic');
                    end
                    if std(tmp)<tst
                        tst = std(tmp);
                        ind = [fac jatf jzeta jww];
                    end
                end    
            end
        end
        jww
    end
end