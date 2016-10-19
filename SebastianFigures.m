% making images for Sebatsian's MHC I paper

para = [75e3/60 [470 670]/1.33 394];

global pd
pd = [0.01 50];

%load 'Y:\FCS\Qui\MHC1 Messungen fürs Manuskript\Ingos Fits\FP7.mat'
load 'Y:\FCS\Qui\ab 091214_wt MHC ohne Peptid\091215_wt MHC ohne Peptid\091215_wt MHC ohne peptid_FRET_5uW.mat'
[data.y,data.t]=FCSCrossRead(res,[1e-6 5]);
data.y(1,:,:)=[]; data.t(1)=[];
%data.y=data.y(:,:,[1:53 55:62 72:end]);
w0 = 350; 
a0 = 150;
[dc v conc w0 a0 triplet c velo err z] = FCSFit(data,[w0 a0],1,0,para);

pd0 = pd;
clear dv 
close all
part1 = 10;
part2 = 20;
for part=part1:part2
    pd = pd0;
    for k=1:floor(size(data.y,3)/part)
        y = mean(data.y(:,:,(k-1)*part+1:k*part),3);
        GaussFcs([w0 a0],para(1),para(2:3),para(4),1e6*data.t,y(:,1:2),(y(:,3)+y(:,4))/2,1,[0 triplet; 1 triplet],[],1./std(data.y(:,1:3,:),0,3));
        dv{part-part1+1}(k) = 1e-8/pd(1);
    end
end
for part=1:part2-part1+1
    if part==1
        tm(part) = mean(res.time);
    else
        tm(part) = mean(sum(reshape(res.time(1:(part1+part-1)*floor(size(data.y,3)/(part1+part-1))),(part1+part-1),floor(size(data.y,3)/(part1+part-1)))));
    end
end
for j=1:part2-part1+1 dcm(j)=mean(dv{j}); dcs(j)=std(dv{j}); rm(j)=mean(Diff2Rad(dv{j},21)); rs(j)=std(Diff2Rad(dv{j},21)); end
h = errorbar(tm,rm,rs); hold on; plot(tm,rm,'ob'); hold off
tmp = get(h,'Children');
tmp = get(tmp(2),'YData');
tmp = tmp(1:9:end);
for part = 1:part2-part1+1
    text(tm(part),tmp(part),int2str(length(dv{part})),'VerticalAlignment','Bottom','HorizontalAlignment','Center')
end
xlabel('measurement time (s)');
ylabel('hydrodynamic radius (Angstrom)')

return

para = [75e3/60 [470 670]/1.33 394];

global pd
pd = [0.01 50];

load 'Y:\FCS\Qui\MHC1 Messungen fürs Manuskript\Ingos Fits\FL9.mat'
[data.y,data.t]=FCSCrossRead(res);
w0 = 350; 
a0 = 150;
[dc v conc w0 a0 triplet c velo err z] = FCSFit(data,[w0 a0],1,0,para);

pd0 = pd;
clear dv 
close all
max_part = min([size(data.y,3),10]);
for part=1:max_part
    pd = pd0;
    for k=1:floor(size(data.y,3)/part)
        y = mean(data.y(:,:,(k-1)*part+1:k*part),3);
        GaussFcs([w0 a0],para(1),para(2:3),para(4),1e6*data.t,y(:,1:2),(y(:,3)+y(:,4))/2,1,[0 triplet; 1 triplet],[],1./std(data.y(:,1:3,:),0,3));
        dv{part}(k) = 1e-8/pd(1);
    end
end
for part=1:max_part
    if part==1
        tm(part) = mean(res.time);
    else
        tm(part) = mean(sum(reshape(res.time(1:part*floor(size(data.y,3)/part)),part,floor(size(data.y,3)/part))));
    end
end
for j=1:length(dv) dcm(j)=mean(dv{j}); dcs(j)=std(dv{j}); rm(j)=mean(Diff2Rad(dv{j},21)); rs(j)=std(Diff2Rad(dv{j},21)); end
h = errorbar(tm,rm,rs); hold on; plot(tm,rm,'ob'); hold off
tmp = get(h,'Children');
tmp = get(tmp(2),'YData');
tmp = tmp(1:9:end);
for part = 1:10
    text(tm(part),tmp(part),int2str(length(dv{part})),'VerticalAlignment','Bottom','HorizontalAlignment','Center')
end
xlabel('measurement time (s)');
ylabel('hydrodynamic radius (Angstrom)')


return

para = [75e3/60 [470 670]/1.33 394];

global pd
pd = [0.01 50];

load 'Y:\FCS\Qui\MHC1 Messungen fürs Manuskript\091103_Y85C_ohne FP7\091103_Y85C ohne FP7_FRET'
[data.y,data.t]=FCSCrossRead(res);
data.y(:,:,[7 9]) = [];

w0 = 350; 
a0 = 150;
[dc v conc w0 a0 triplet c velo err z] = FCSFit(data,[w0 a0],1,0,para);

pd0 = pd;
clear dv dcm dcs rm rs tm
close all
max_part = min([size(data.y,3),10]);
for part=1:max_part
    pd = pd0;
    for k=1:floor(size(data.y,3)/part)
        y = mean(data.y(:,:,(k-1)*part+1:k*part),3);
        GaussFcs([w0 a0],para(1),para(2:3),para(4),1e6*data.t,y(:,1:2),(y(:,3)+y(:,4))/2,1,[0 triplet; 1 triplet],[],1./std(data.y(:,1:3,:),0,3));
        dv{part}(k) = 1e-8/pd(1);
    end
end
for part=1:max_part
    if part==1
        tm(part) = mean(res.time);
    else
        tm(part) = mean(sum(reshape(res.time(1:part*floor(size(data.y,3)/part)),part,floor(size(data.y,3)/part))));
    end
end
for j=1:length(dv) dcm(j)=mean(dv{j}); dcs(j)=std(dv{j}); rm(j)=mean(Diff2Rad(dv{j},21)); rs(j)=std(Diff2Rad(dv{j},21)); end
h = errorbar(tm,rm,rs); hold on; plot(tm,rm,'ob'); hold off
tmp = get(h,'Children');
tmp = get(tmp(2),'YData');
tmp = tmp(1:9:end);
for part = 1:max_part
    text(tm(part),tmp(part),int2str(length(dv{part})),'VerticalAlignment','Bottom','HorizontalAlignment','Center')
end
xlabel('measurement time (s)');
ylabel('hydrodynamic radius (Angstrom)')
