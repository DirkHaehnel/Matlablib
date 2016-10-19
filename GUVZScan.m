% GUV-2fFcs

path = 'D:\Doc\FCS\2Focus\GUV\2006-09-07GUV_zscan\';
%path = 'D:\Joerg\Doc\Fcs\2Focus\GUV\2006-09-05 GUV_z-scan\';
%path = 'D:\Joerg\Doc\Fcs\2Focus\GUV\2006-09-06_GUV_leistungsabhängig\';
names = dir([path '*c.mat'])
close all

global pd

if 1
    clear dc v conc w0 a0 triplet c
    for j=1:length(names)
        for cnt = 1:1
            %pd = [0.1 5e4];
            pd = 0.1;
            [dc(j,cnt) v(j,cnt) conc(j,cnt) w0(:,j,cnt)] = FCSFit([path names(j).name],[400 400 i]);
            %[dc(j,cnt) v(j,cnt) conc(j,cnt) w0(:,j,cnt) triplet(:,j)] = FCSFit([path names(j).name],[400 400 i]);
            %[dc(:,j,cnt) v(j,cnt) conc(j,cnt) w0(:,j,cnt)] = FCSFit([path names(j).name],[400 400 i],0);
        end
        %eval(['print -dpng -r300 ''' names(j).name(1:end-4) ''''])
    end
    %close; p = Simplex('LaserBeamFit',[300 0],[],[],[],[],640/1.33,(zb(2:8)-40.5)*1e3,w0(2:8),[],1);
    %         close; p = Simplex('LaserBeamFit',[300 0],[],[],[],[],640/1.33,(zc-40.5)*1e3,w0,[],1);
    %         for k=1:5
    %             p = Simplex('LaserBeamFit',p,[],[],[],[],640/1.33,(zc-40.5)*1e3,w0,[],1);
    %         end
    %         [err , c] = LaserBeamFit(p,640/1.33,(zc-40.5)*1e3,w0,[],1);
end

if 0
    clear p1 p2
    for j=1:length(names)
        pd = [0.1 10e3];
        [y,t] = FCSCrossRead([path names(j).name]);
        pb1(:,j) = Simplex('GaussFCS',400,[],[],[],[],100e3/60,[640 670]/1.33,[],1e6*t,sum(y(:,1,isfinite(sum(sum(y(:,1,:),1),2))),3),[],1,[[0.1 0];[0.1 inf]]);
        pb2(:,j) = Simplex('GaussFCS',400,[],[],[],[],100e3/60,[640 670]/1.33,[],1e6*t,sum(y(:,2,isfinite(sum(sum(y(:,1,:),1),2))),3),[],1,[[0.1 0];[0.1 inf]]);
    end
end