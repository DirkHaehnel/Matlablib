% program DsomFcs tries to calculate reduced-volume FCS by using triplet
% state dynamics

t_trip = 10e3; % 10 mysec
timewin = 1; % time window of data processing in s
timeframe = 3e10; % time-frame to read photons from file in ns
smooth = 1; % smoothed autocorrelation
close all

cd D:\Joerg\Doc\Sykora\DSOM
name1 = 'D:\Joerg\Doc\Sykora\DSOM\Al647_f050_chb1.pt3';
name2 = 'D:\Joerg\Doc\Sykora\DSOM\Al647_f050_chb2.pt3';

% [hofcs,head] = pt3pro(name1, name2, 'hofcs');
% save AlexaHoFCS hofcs head

% [res, head] = pt3pro(name1, name2, 'xfcs');
% save AlexaFCS res head

load AlexaFCS 

sync = 1e9/head.CntRate0;

Nsub = 10;
Ncasc = ceil(log2(3/(t_trip*1e-9*Nsub)+1));,

triptime = (0:t_trip/sync)';
autotime = 1e-9*t_trip*cumsum(reshape(repmat(2.^[0 0:Ncasc-1],Nsub,1),(Ncasc+1)*Nsub,1));

p_trip = [3e5 0];
len = length(p_trip);
auto = zeros(length(autotime),len,len);
para = exp(-sync*1e-9*triptime*p_trip);
% filter generation
weight = interp1(res.autotime,mean(res.auto,2),triptime*sync*1e-9,'cubic'); 
for j=1:len
    para(:,j) = para(:,j)/sum(para(:,j));
end % normalization
master = (para'*diag(1./weight)*para)\(diag(1./weight)*para)';
master = master';

head = pt3Read_t(name1);
t_Counts = head.NCounts;
head = pt3Read_t(name2);
t_Counts = t_Counts + head.NCounts;

rate    = [];
time    = [];
tst1    = 1;
tst2    = 1;
cnt     = 0;
cnt1    = 0;
cnt2    = 0;

max_time  = 0;

overcount1 = 0;
overcount2 = 0;

t_trip = t_trip/sync;

while ((tst1>0)&(tst2>0))

    max_time = max_time + timeframe;

    [head, tmpy1, tmpx1, tst1, overcount1] = pt3Read_t(name1, cnt1+1, overcount1, max_time);
    [head, tmpy2, tmpx2, tst2, overcount2] = pt3Read_t(name2, cnt2+1, overcount2, max_time);

    cnt1 = cnt1 + tst1;
    cnt2 = cnt2 + tst2;
    
    if ((tst1>0)&(tst2>0))

        tmpy = ([tmpy1+tmpx1; tmpy2+tmpx2]/sync);  % in timeunits

        ind = [zeros(length(tmpy1),1); ones(length(tmpy2),1)];

        [tmpy, ord] = sort(tmpy);
        ind = ind(ord);
        tmpy = tmpy-tmpy(1);
        ind1 = zeros(length(ind),len);
        ind2 = ind1;
        t0 = (1:length(ind))';
        
        for k=1:1
           tmp = [round(tmpy(k+1:end)-tmpy(1:end-k)); ones(k,1)*2*t_trip];
           t = t0(tmp<t_trip);
           for j=1:len
              ind1(t,j) = ind1(t,j) + double(ind(t)==0 & ind(t+k)==1).*master(tmp(t)+1,j);
              ind2(t,j) = ind2(t,j) + double(ind(t)==1 & ind(t+k)==0).*master(tmp(t)+1,j);              
           end
        end

        cnt = cnt+length(tmpy);
        dt = (tmpy(end))*sync*1e-9;
        time = [time dt];
        
        tmp = ind1(:,1)==0 & ind2(:,1)==0;
        tmpy(tmp) = [];
        ind1(tmp,:) = [];
        ind2(tmp,:) = [];

        for j1=1:size(para,2)
            for j2=1:size(para,2)
                auto(:,j1,j2) = auto(:,j1,j2) + tttr2fcs(round(tmpy/t_trip), ind1(:,j1), ind2(:,j2), Ncasc, Nsub, smooth);
                auto(:,j1,j2) = auto(:,j1,j2) + tttr2fcs(round(tmpy/t_trip), ind2(:,j1), ind1(:,j2), Ncasc, Nsub, smooth);
            end
        end

        subplot('position', [0.875 0.2 0.05 0.6]);
        bar(0,(cnt1+cnt2)/t_Counts);
        axis([-0.4 0.4 0 1]);
        set(gca,'xtick',[],'ytick',[]);
        subplot('position', [0.1 0.1 0.7 0.8]);
        tmp = auto(2:end,1,1);
        for j=2:len 
            tmp = [tmp auto(2:end,j,j)];
        end
        FCSScalePlot(autotime(2:end),tmp/sum(time),1);
        drawnow;
    end

end
clf

tmp = auto(2:end,1,1);
for j=2:len
    tmp = [tmp auto(2:end,j,j)];
end
FCSScalePlot(autotime(2:end),tmp/sum(time),1);

save AlexaDSOMFcs auto autotime

return

p1=Simplex('Rigler',[1e-2 1e-3 1e4],[],[],[],[],autotime(2:end),auto(2:end,1,1));
for k=1:2
    p1=Simplex('Rigler',p1,[],[],[],[],autotime(2:end),auto(2:end,1,1));
end
p2=Simplex('Rigler',p1,[],[],[],[],autotime(2:end),auto(2:end,2,2),[],1)
for k=1:2
    p2=Simplex('Rigler',p2,[],[],[],[],autotime(2:end),auto(2:end,2,2),[],1);
end

[err c1]=Rigler(p1,autotime(2:end),auto(2:end,1,1),[],1);
[err c2]=Rigler(p2,autotime(2:end),auto(2:end,2,2),[],1);
[err c1 z1]=Rigler(p1,c1,autotime(2:end));
[err c2 z2]=Rigler(p2,c2,autotime(2:end));
semilogx(autotime(2:end),(auto(2:end,1,1)-c1(1))/c1(2)/1.09,'or',autotime(2:end),(auto(2:end,2,2)-c2(1))/c2(2),'db',autotime(2:end),(z1*c1-c1(1))/c1(2)/1.09,'r',autotime(2:end),(z2*c2-c2(1))/c2(2),'b','markersize',4)
axis([2e-5 3 0 1.6])
xlabel('time [s]'); ylabel('autocorrelation [a.u.]')
legend({'standard ACF','DSOM-ACF'})

