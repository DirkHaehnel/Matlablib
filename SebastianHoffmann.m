% program Seastian Hoffmann for modelling particle interaction FCS

if 0 % calculate MDF alias PSF
    clear all
    close all
    rhofield =  [0 1];
    focpos = 45;
    zfield = focpos + [0 2];
    NA = 0.9;
    n0 = 1.497; % toluol
    n = n0;
    n1 = n0;
    d0 = [];
    d = 0;
    d1 = [];
    mag = 40;
    fd = 164.5e3/mag;
    lamex = 0.543;
    over = 5e3;
    ring = [];
    av = 70/2;
    lamem = 0.570;
    zpin = 0;

    exc = GaussExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos,[],lamex/0.005);
    mdf = GaussExc2MDF(exc, NA, n0, n, n1, focpos, lamem, mag, av, zpin);
    [modres, autotime] = FCS(exc.rho,exc.z,mdf.volx,mdf.voly);

    u = squeeze(mdf.volx(:,:,1) + mdf.voly(:,:,1));

    save SebastianHoffmannModel exc mdf u modres autotime
end



if 1
    load SebastianHoffmannModel
    
    auto = FCSInteractionDH([0.1 1],10.^(-6:0.1:0),exc.rho,exc.z,u);
    
    return
    
    name = 'c:\Joerg\Doc\Fcs\FCSInteraction\PS75k_20wt.fcs';
    name = 'c:\Joerg\Doc\Fcs\FCSInteraction\PS75k_33wt.fcs';
    fin = fopen(name,'r');
    tst = true;
    while tst
        if strcmp(fgetl(fin),'##XYPOINTS=(XY..XY)')
            tst = false;
        end
    end
    x = fscanf(fin,'%f, %f\n', inf);
    fclose(fin)
    t = x(1:2:end);
    x = x(2:2:end);
    x = x(t>1e-3);
    t = t(t>1e-3);

    close; p=Simplex('FCSInteractionFun',[0.5 3],[0 0],[],[],[],t*1e3,x-1,autotime,modres,exc.rho,exc.z-exc.focpos(1),mdf.volx(:,:,1)+mdf.voly(:,:,1));
end