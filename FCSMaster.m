cd D:\Joerg\Doc\Gregor\DiffStandards\040216

if 1
    clear all
%    load FCSMODEL
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
    overv = 5.3e3;
    focpos = 10;
    lamem = 0.67;
    mag = tubelens/wd; %60;
    av = 100/2; %pinhole radius!!!
    zpinv = 0:500; 
    atfv = 0;
    kappa = 0; 
    lt = [];
    satv = 0:0.1:10;
    pow = [1 3 10 30 100 300]/300;
    %pow = [0.2000    0.6670    2.0000    6.6700   20.0000   66.7000  200.0000]/200;
    
%     jo0 = jo;
%     jz0 = jz;
%     ja0 = ja;
%     js0 = js; 
    jo0 = 1;
    jz0 = 1;
    ja0 = 1;
    js0 = 1; 

    for jo=jo0:length(overv)
        over = [wd*NA inf overv(jo)];        
        if jo==jo0 jz1 = jz0; else jz1 = 1; end
        for jz=jz1:length(zpinv)
            zpin = zpinv(jz);
            if (jo==jo0 & jz==jz0) ja1 = ja0; else ja1 = 1; end            
            for ja=ja1:length(atfv)   
                atf = [1.51 atfv(ja)];
                if (jo==jo0 & jz==jz0 & ja==ja0) js1 = js0; else js1 = 1; end                            
                for js=js1:length(satv)
                    sat = satv(js);
                    [jo jz ja js]
                    [modres(:,jo,jz,ja,js,:), autotime, veff(jo,jz,ja,js,:) exc(:,jo,jz,ja,js,:)] = ...
                        NewFCS(NA, n0, n, n1, d0, d, d1, lamex, lamem, mag, focpos, av, over, zpin, kappa, atf, lt, sat*pow);
                end
                save FCSMODEL
            end
        end
    end
end

if 1
    clear all
    close all
    load alexa647
    load FCSMODEL modres autotime
    t = res1.autotime;
    for j=1:6
        eval(['x(:,' num2str(j) ') = mean(res' num2str(7-j) '.auto'')'';']);
    end
    tt = 10:length(t);
    jo0 = 1;
    jz0 = 1;
    ja0 = 1;
    js0 = 1; 
    for jo=jo0:size(modres,2)
        if jo==jo0 jz1 = jz0; else jz1 = 1; end
        for jz=jz1:size(modres,3)
            if (jo==jo0 & jz==jz0) ja1 = ja0; else ja1 = 1; end            
            for ja=ja1:size(modres,4)
                if (jo==jo0 & jz==jz0 & ja==ja0) js1 = js0; else js1 = 1; end                            
                for js=js1:size(modres,5)
                    p(:,jo,jz,ja,js) = simplex('alexafit',[5e-7 1e6 1e4 100],zeros(1,4),[],[],[],pow,t(tt),x(tt,:),autotime,squeeze(modres(:,jo,jz,ja,js,1:6)));
                    % p(jo,jz,ja,js) = simplex('affineexpfit',[1e-5 1],[0 0],[],[],[],t(tt),x(tt,:),autotime,squeeze(modres(:,jo,jz,ja,js,:)));
                    [err(jo,jz,ja,js), bla, z(:,jo,jz,ja,js,:)] = alexafit(p(:,jo,jz,ja,js),pow,t(tt),x(tt,:),autotime,squeeze(modres(:,jo,jz,ja,js,1:6)));
                end
            end
        end
    end
    
    tst = err(1,1,1,1);
    ind = [1 1 1 1];
    jo0 = 1;
    jz0 = 1;
    ja0 = 1;
    js0 = 1; 
    for jo=jo0:size(modres,2)
        if jo==jo0 jz1 = jz0; else jz1 = 1; end
        for jz=jz1:size(modres,3)
            if (jo==jo0 & jz==jz0) ja1 = ja0; else ja1 = 1; end            
            for ja=ja1:size(modres,4)
                if (jo==jo0 & jz==jz0 & ja==ja0) js1 = js0; else js1 = 1; end                            
                for js=js1:size(modres,5)
                    if err(jo,jz,ja,js)<tst
                        tst = err(jo,jz,ja,js);
                        ind = [jo,jz,ja,js];
                    end
                end
            end
        end
    end
    save fitresults z err p ind
end
    

if 0
    ind = repmat([1 1 1 1],6,1);
    jo0 = 1;
    jz0 = 1;
    ja0 = 1;
    js0 = 1; 
    for j=1:6
        eval(['err = err' num2str(j) '(1,1,1,1);'])
        tst(j) = err;
    end
    for jo=8:8%jo0:size(modres,2)
        if jo==jo0 jz1 = jz0; else jz1 = 1; end
        for jz=2:2%jz1:size(modres,3)
            if (jo==jo0 & jz==jz0) ja1 = ja0; else ja1 = 1; end            
            for ja=9:9%ja1:size(modres,4)
                if (jo==jo0 & jz==jz0 & ja==ja0) js1 = js0; else js1 = 1; end                            
                for js=js1:size(modres,5)
                    err(1) = err1(jo,jz,ja,js);
                    err(2) = err2(jo,jz,ja,js);
                    err(3) = err3(jo,jz,ja,js);
                    err(4) = err4(jo,jz,ja,js);
                    err(5) = err5(jo,jz,ja,js);
                    err(6) = err6(jo,jz,ja,js);
                    p(1) = p1(jo,jz,ja,js);
                    p(2) = p2(jo,jz,ja,js);
                    p(3) = p3(jo,jz,ja,js);
                    p(4) = p4(jo,jz,ja,js);
                    p(5) = p5(jo,jz,ja,js);
                    p(6) = p6(jo,jz,ja,js);
                    p = 5e-13./p*1e7;
                    for j=1:6
                        if err(j)<tst(j) & p(j)>1.3 & p(j)<1.5
                            tst(j) = err(j);
                            ind(j,:) = [jo,jz,ja,js];
                        end
                    end
                end
            end
        end
    end
end


if 0
    natconst
    T0 = 273.15 + 26;
    load('D:\Joerg\Doc\Gregor\DiffStandards\viscosity.mat');
    eta = interp1(T,visc(1,:),T0,'cubic')*1e-3;
    rad = 36e-9/2;
    1e4*BoltzmannConstant*T0/6/pi/eta/rad
end