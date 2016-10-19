close all
clear all

av = 100e3/60;
lam = [640 670]/1.33;
dist = 403;

if 0
    names = dir([dirname '*.t3r']);
    for j=1:length(names)
        if exist(['c:\Joerg\Doc\Fcs\2Focus\Recoverin\Vesicle\' names(j).name(1:end-4) '.mat'])==0
            [res head] = TwoFocus23FCS([dirname names(j).name]);
            save([dirname names(j).name(1:end-4)],'res','head');
        end
    end
end

if 0
    dirnamelist = {'c:\Joerg\Doc\Fcs\2Focus\Recoverin\071107_10mgmlDOPC',...
        'c:\Joerg\Doc\Fcs\2Focus\Recoverin\071112_5mgmlDOPC',...
        'c:\Joerg\Doc\Fcs\2Focus\Recoverin\20mgmlROSMix',...
        'c:\Joerg\Doc\Fcs\2Focus\Recoverin\m-2_20mgmlDOPC',...
        'c:\Joerg\Doc\Fcs\2Focus\Recoverin\m-3_20mgmlDOPC',...
        'c:\Joerg\Doc\Fcs\2Focus\Recoverin\ohneChelator',...
        'c:\Joerg\Doc\Fcs\2Focus\Recoverin\wtmRec_20mgmlDOPC'};
    for jdir=1:length(dirnamelist)
        dirname = dirnamelist{jdir};
        clear fitres0;
        names = dir([dirname '\*3photon.mat']);
        for j=1:length(names)
            load([dirname '\' names(j).name(1:end-4)]);

            clear data

            ind = res.autotime>1e-5 & res.autotime<1;
            res.auto = res.auto(ind,:,:,:);
            res.auto2 = res.auto2(ind,:,:,:);
            res.auto3 = res.auto3(ind,:,:,:,:);
            res.autotime = res.autotime(ind);
            
            ind = abs(mean(cov(squeeze(res.auto(:,1,3,1:end))))-mean(mean(cov(squeeze(res.auto(:,1,3,1:end))))))>std(mean(cov(squeeze(res.auto(:,1,3,1:end)))));
            res.auto(:,:,:,ind) = [];
            res.auto2(:,:,:,ind) = [];
            res.auto3(:,:,:,:,ind) = [];
            res.rate(ind,:) = [];


            data.y(:,1,:) = sqrt(res.auto(:,1,3,:).*res.auto(:,3,1,:));
            data.y(:,2,:) = sqrt(res.auto(:,2,4,:).*res.auto(:,4,2,:));
            data.y(:,3,:) = sqrt(res.auto(:,1,4,:).*res.auto(:,3,2,:));
            data.y(:,4,:) = sqrt(res.auto(:,2,3,:).*res.auto(:,4,1,:));
            data.t = res.autotime;
            

            global pd
            %             pd = [1e-8./[1e-6 1e-7] 10 20];
            %             bounds = [1e-8./2e-6 0 0 0; inf 1e-8./2e-8 100 100];
            pd = [1e-8./[1e-6 1e-7] 10];
            bounds = [1e-8./2e-6 0 0; inf 1e-8./2e-8 100];
            [dc v conc w a triplet c velo err z] = FcsFit(data,[350 150],length(pd)-2,[],[av lam dist],bounds);
            [ind, ind] = sort(dc);
            pd(1:2) = pd(ind);
            [dc v conc w a triplet c velo err z] = FcsFit(data,[w a],length(pd)-2,[],[av lam dist],bounds);
            fitres0(j).dc = dc;
            fitres0(j).w1 = w;
            fitres0(j).a1 = a;
            fitres0(j).c1 = c;
            fitres0(j).triplet1 = triplet;
            fitres0(j).z1 = z;
            fitres0(j).err1 = err;

            [zd zd] = TripleFCS(res); zd = sqrt(prod(sum(zd,3),2));
            
            close all
            ind = res.autotime>1e-5;
            %             bounds = [pd(1:2)' 0 0; pd(1:2)' 100 100];
            %             pd = [pd(1:2)' 10 20];
            bounds = [pd(1:2)' 0; pd(1:2)' 100];
            pd = [pd(1:2)' 20];
            para = Simplex('Gauss3Fcs',[w a],[],[],[],[],av,lam,res.autotime(ind)*1e6,zd(ind,:),1,bounds);
            para = Simplex('Gauss3Fcs',para,[],[],[],[],av,lam,res.autotime(ind)*1e6,zd(ind,:),1,bounds);
            [err, c2, z2] = Gauss3Fcs(para,av,lam,res.autotime(ind)*1e6,zd(ind,:),1,bounds);
            fitres0(j).c2 = c2;
            fitres0(j).z2 = z2;
            fitres0(j).w2 = para(1);
            fitres0(j).a2 = para(2);
            fitres0(j).err2 = err;

            eval(['save ' dirname '\IrisVesicleFitres fitres0'])
        end
    end
    % for j=1:length(fitres0) [fitres0(j).dc ind]=sort(fitres0(j).dc); tmp=fitres0(j).c1(2:3,:); tmp=tmp(ind,:); fitres0(j).c1(2:3,:)=tmp; end
    % caconc = [65e-9 180e-9 380e-9 830e-9 2.8e-6 2.1e-6 9.5e-6 35e-6 2.1e-3 20e-9 120e-9 290e-9 545e-9 1.35e-6 18e-6 3e-6 21e-6 513e-6];
    % save IrisVesicleFitres fitres0 caconc -append
end

if 0 % raw fit inspection
    dirname = 'd:\071107_10mgmlDOPC';
    %     dirname = 'c:\Joerg\Doc\Fcs\2Focus\Recoverin\071107_10mgmlDOPC';
    %     dirname = 'c:\Joerg\Doc\Fcs\2Focus\Recoverin\071112_5mgmlDOPC';
    %     dirname = 'c:\Joerg\Doc\Fcs\2Focus\Recoverin\20mgmlROSMix';
    %     dirname = 'c:\Joerg\Doc\Fcs\2Focus\Recoverin\m-2_20mgmlDOPC';
    %     dirname = 'c:\Joerg\Doc\Fcs\2Focus\Recoverin\m-3_20mgmlDOPC';
    %     dirname = 'c:\Joerg\Doc\Fcs\2Focus\Recoverin\ohneChelator';
    %     dirname = 'c:\Joerg\Doc\Fcs\2Focus\Recoverin\wtmRec_20mgmlDOPC';
    names = dir([dirname '\*3photon.mat']);
    clear pnum
    for j=1:length(names)
        tmp=findstr(names(1).name,'_');
        pnum(j) = str2num(names(j).name(tmp(end-1)+2:tmp(end)-1));
    end
    caconc = [20e-9; 65e-9; 120e-9; 180e-9; 290e-9; 380e-9; 545e-9; 830e-9; 1.35e-6; 2.8e-6; 18e-6; 2.1e-6; 3e-6; 9.5e-6; 21e-6; 35e-6; 513e-6; 2.1e-3];    
    
	load c:\Joerg\Doc\Fcs\2Focus\Recoverin\071107_10mgmlDOPC\IrisVesicleFitres.mat
    for j=1:length(fitres0)
        dc(:,j) = fitres0(j).dc;
        c1(:,j) = mean(fitres0(j).c1(2:3,:),2);
        c2(:,j) = fitres0(j).c2(2:3);
        cc(:,j) = c1(:,j).^3./c2(:,j).^2;
    end
    loglog(caconc(pnum),dc,'o-',caconc(pnum),ones(1,length(pnum))*mean(dc(1,end-5:end)),caconc(pnum),ones(1,length(pnum))*mean(dc(2,end-5:end)))
end

if 1 % fit with fixed diffusion coefficients
    % dirname = 'd:\071107_10mgmlDOPC';
    % dirname = 'd:\071112_5mgmlDOPC';
    %     dirname = 'c:\Joerg\Doc\Fcs\2Focus\Recoverin\20mgmlROSMix';
    dirname = 'd:\m-2_20mgmlDOPC';
    % dirname = 'd:\m-3_20mgmlDOPC';
    %     dirname = 'c:\Joerg\Doc\Fcs\2Focus\Recoverin\ohneChelator';
    %     dirname = 'c:\Joerg\Doc\Fcs\2Focus\Recoverin\wtmRec_20mgmlDOPC';
    pd0 = 1e-8./1e-6;
    names = dir([dirname '\*250C*3photon.mat']);
    % clear fitres
    for j=1:length(names)
        load([dirname '\' names(j).name(1:end-4)]);

        clear data

        ind = res.autotime>1e-5 & res.autotime<10;
        res.auto = res.auto(ind,:,:,:);
        res.auto2 = res.auto2(ind,:,:,:);
        res.auto3 = res.auto3(ind,:,:,:,:);
        res.autotime = res.autotime(ind);

        ind = abs(mean(cov(squeeze(res.auto(:,1,3,1:end))))-mean(mean(cov(squeeze(res.auto(:,1,3,1:end))))))>std(mean(cov(squeeze(res.auto(:,1,3,1:end)))));
        res.auto(:,:,:,ind) = [];
        res.auto2(:,:,:,ind) = [];
        res.auto3(:,:,:,:,ind) = [];
        res.rate(ind,:) = [];

        data.y(:,1,:) = sqrt(res.auto(:,1,3,:).*res.auto(:,3,1,:));
        data.y(:,2,:) = sqrt(res.auto(:,2,4,:).*res.auto(:,4,2,:));
        data.y(:,3,:) = sqrt(res.auto(:,1,4,:).*res.auto(:,3,2,:));
        data.y(:,4,:) = sqrt(res.auto(:,2,3,:).*res.auto(:,4,1,:));
        data.t = res.autotime;

        global pd
        pd = [pd0 0.05 20];
        bounds = [pd0 0 0; pd0 inf 100];
        [dc v conc w a triplet c velo err z] = FcsFit(data,[420 190],length(pd)-2,[],[av lam dist],bounds);
        fitres(j).dc = dc;
        fitres(j).w1 = w;
        fitres(j).a1 = a;
        fitres(j).c1 = c;
        fitres(j).triplet1 = triplet;
        fitres(j).z1 = z;
        fitres(j).err1 = err;

        [zd zd] = TripleFCS(res); zd = sum(sum(zd,3),2);

        close all
        bounds = [pd0 pd(2) 0 0; pd0 pd(2) 100 100];
        pd = [pd0 pd(2) 10 10];
        para = Simplex('Gauss3Fcs',[w a],[],[],[],[],av,lam,res.autotime*1e6,zd,1,bounds);
        para = Simplex('Gauss3Fcs',para,[],[],[],[],av,lam,res.autotime*1e6,zd,1,bounds);
        [err, c2, z2] = Gauss3Fcs(para,av,lam,res.autotime*1e6,zd,1,bounds);
        fitres(j).c2 = c2;
        fitres(j).z2 = z2;
        fitres(j).w2 = para(1);
        fitres(j).a2 = para(2);
        fitres(j).err2 = err;

        save IrisVesicleFitres250 fitres
    end
end

if 0
    dirname = 'c:\Joerg\Doc\Fcs\2Focus\Recoverin\071107_10mgmlDOPC';
    %     dirname = 'c:\Joerg\Doc\Fcs\2Focus\Recoverin\071112_5mgmlDOPC';
    %     dirname = 'c:\Joerg\Doc\Fcs\2Focus\Recoverin\20mgmlROSMix';
    %     dirname = 'c:\Joerg\Doc\Fcs\2Focus\Recoverin\m-2_20mgmlDOPC';
    %     dirname = 'c:\Joerg\Doc\Fcs\2Focus\Recoverin\m-3_20mgmlDOPC';
    %     dirname = 'c:\Joerg\Doc\Fcs\2Focus\Recoverin\ohneChelator';
    %     dirname = 'c:\Joerg\Doc\Fcs\2Focus\Recoverin\wtmRec_20mgmlDOPC';
    names = dir([dirname '\*3photon.mat']);
    clear pnum
    for j=1:length(names)
        tmp=findstr(names(1).name,'_');
        pnum(j) = str2num(names(j).name(tmp(end-1)+2:tmp(end)-1));
    end
    caconc = [20e-9; 65e-9; 120e-9; 180e-9; 290e-9; 380e-9; 545e-9; 830e-9; 1.35e-6; 2.8e-6; 18e-6; 2.1e-6; 3e-6; 9.5e-6; 21e-6; 35e-6; 513e-6; 2.1e-3];    
    
	load([dirname '\IrisVesicleFitres' ]);
    clear dc c1 c2 cc c10 c20 cc0
    for j=1:length(fitres)
        dc(:,j) = fitres(j).dc;
%         c10(:,j) = mean(fitres0(j).c1(end-1:end,:),2);
%         c20(:,j) = fitres0(j).c2(end-1:end);
%         cc0(:,j) = c10(:,j).^3./c20(:,j).^2;
        c1(:,j) = mean(fitres(j).c1(end-1:end,:),2);
        c2(:,j) = fitres(j).c2(end-1:end);
        cc(:,j) = c1(:,j).^3./c2(:,j).^2;
    end
    [conc,ind] = sort(caconc(pnum));
    close; phill = Simplex('Hill',[1e-6 3], [],[],[],[],  conc, cc(2,ind)./sum(cc(:,ind)));
    xlabel('Ca conc. (M)')
    ylabel('bound part')

    for jdir=1:length(dirnamelist)
        dirname = dirnamelist{jdir};
        load([dirname '\IrisVesicleFitres' ]);
        for j=1:length(fitres)
            c1(:,j) = mean(fitres(j).c1(2:3,:),2);
            c2(:,j) = fitres(j).c2(2:3);
            cc(:,j) = c1(:,j).^3./c2(:,j).^2;
        end
        [ind,ind]=sort(caconc);
        close; phill = Simplex('Hill',[1e-6 3], [],[],[],[],  caconc(ind), cc(1,ind)./sum(cc(:,ind)));
        xlabel('Ca conc. (M)')
        ylabel('bound part')
    end
end

return
% old routines

path = 'D:\Doc\Iris\NewVesicles\';
name = 'wtmRecX647_nachZentrifuge_250C.mat';

global pd
pd = [0.01 1e5];
FCSFit([path name]);
pd0 = pd(1);

pd = [1e-4; pd];

files = dir([path '*_P*.mat'])
for j=1:length(files)
    pos = findstr(files(j).name,'_P')+2;
    t(j) = str2num(files(j).name(pos:pos+1));
    [dc v conc w0 a0 triplet c err] = FCSFit([path files(j).name],[],1,[],[],[[0 pd0 0];inf pd0 inf]);
    dc1(j) = dc(1); dc2(j) = dc(2); % dc1 = vector of slow diffusion coefficient
    c1(:,j) = c(:,1); c2(:,j) = c(:,2); % c1, c2 = amplitudes of different components in ACF 1 and ACF 2
end
    
plot(t,dc1)
xlabel('number of measurement')
ylabel('slow diff. coef.')

figure
plot(t,mean([c1(2,:)./sum(c1(2:3,:));c2(2,:)./sum(c2(2:3,:))]))
xlabel('number of measurement')
ylabel('\chi_{slow}')


return

caconc = [20e-9; 65e-9; 120e-9; 180e-9; 290e-9; 380e-9; 545e-9; 830e-9; 1.35e-6; 2.8e-6; 18e-6; 2.1e-6; 3e-6; 9.5e-6; 21e-6; 35e-6; 513e-6; 2.1e-3];
