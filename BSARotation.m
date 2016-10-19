function BSARotation(name)
%close all

% name = 'X:\MTusers\Christina\rotational diffusion\090909_ovalbumin_21C_DICpr_out.ht3';
% name = 'X:\MTusers\Christina\rotational diffusion\090912_aldolase_in_aldol200nM_plusPBS_208C_DICpr_out.ht3';
% name = 'X:\MTusers\Christina\rotational diffusion\091013_HSA_208C_DICpr_out.ht3';

%name = 'm:\MTusers\Christoph\2010-05-10 Amylase Cy5\Amylase Cy5 bis RotDiff 200uW.ht3';

head = ht3v2read_head(name);

if ((isempty(head)) && not(isempty(strfind(name, '.ht3'))))
    % no header => old ht3 File
    fprintf(1,'\n\n      Warning: No File Header. Assuming old ht3 File\n');
    head.Ident = 'HydraHarp';
    head.Resolution = 2e-3;     % in nanoseconds
    head.SyncRate = 80e6;       % 80 MHz sync Rate
    head.DataStart = 0;
    tmp = dir(name);
    head.Records = tmp.bytes/4;
end;

photons = 1e6;
microresolution = head.Resolution;
syncperiod = (1 / head.SyncRate) * 1e9;
automaxtime = 3e3;
frame = round([1 0.5 5]/microresolution);

NCounts = head.Records;

[sync, micro, chan] = ht3v2read_old(name, [1, photons]);
[tcspcraw bin] = mHist(micro(chan==0),0:2^15-1);
tcspcraw(:,2) = mHist(micro(chan==1),0:2^15-1);

[tmp tmp] = max(tcspcraw(:,1));
pos(1) = tmp(1);
ind = microresolution*(bin-bin(pos(1)))>-0.5 & microresolution*(bin-bin(pos(1)))<5;
[tmp tmp] = max(tcspcraw(:,1).*(~ind'));
pos(2) = tmp(1);
pos = sort(pos);
timeshift = zeros(length(bin),2);
weight = zeros(length(bin),2);
timeshift(bin<bin(pos(1))-frame(1),1) = microresolution*bin(pos(2)) - syncperiod;
timeshift(bin>bin(pos(1))+frame(2) & bin<bin(pos(1))+frame(3),1) = microresolution*bin(pos(1));
timeshift(bin>bin(pos(2))+frame(2) & bin<bin(pos(2))+frame(3),1) = microresolution*bin(pos(2));
weight(bin<bin(pos(1))-frame(1),1) = 2;
weight(bin>bin(pos(1))+frame(2) & bin<bin(pos(1))+frame(3),1) = 1;
weight(bin>bin(pos(2))+frame(2) & bin<bin(pos(2))+frame(3),1) = 2;

[tmp tmp] = max(tcspcraw(:,2));
pos(1) = tmp(1);
ind = microresolution*(bin-bin(pos(1)))>-0.5 & microresolution*(bin-bin(pos(1)))<5;
[tmp tmp] = max(tcspcraw(:,2).*(~ind'));
pos(2) = tmp(1);
pos = sort(pos);
timeshift(bin<bin(pos(1))-frame(1),2) = microresolution*bin(pos(2)) - syncperiod;
timeshift(bin>bin(pos(1))+frame(2) & bin<bin(pos(1))+frame(3),2) = microresolution*bin(pos(1));
timeshift(bin>bin(pos(2))+frame(2) & bin<bin(pos(2))+frame(3),2) = microresolution*bin(pos(2));
weight(bin<bin(pos(1))-frame(1),2) = 2;
weight(bin>bin(pos(1))+frame(2) & bin<bin(pos(1))+frame(3),2) = 1;
weight(bin>bin(pos(2))+frame(2) & bin<bin(pos(2))+frame(3),2) = 2;

if (length(unique(timeshift(~timeshift(:,2)==0,2))) ~= length(unique(timeshift(~timeshift(:,1)==0,1)))) 
    tt = min(timeshift, 2);
    ind = tt < 0;
    timeshift(ind,1) = 0;
    timeshift(ind,2) = 0;
    weight(ind,1) = 0;
    weight(ind,2) = 0;
end

autotime = unique(timeshift(~timeshift(:,2)==0,2))*ones(1,length(unique(timeshift(~timeshift(:,2)==0,2))))-...
    ones(length(unique(timeshift(~timeshift(:,1)==0,1))),1)*unique(timeshift(~timeshift(:,1)==0,1))';
autotime = [autotime(:); -autotime(:)];
len = length(autotime);
for j=1:ceil(automaxtime/syncperiod)+1;
    autotime = [autotime; autotime(end-len+1:end)+syncperiod];
end
% autotime(autotime<0 | autotime>automaxtime) = [];
% autotime = sort(unique([0; autotime; max(autotime)+frame(3)*microresolution]));
autotime(autotime>automaxtime) = [];
autotime = unique([min(autotime)-frame(3)*microresolution; autotime; max(autotime)+frame(3)*microresolution]);

auto = 0*autotime;
auto = repmat(auto,[1 2 2 2 ceil(NCounts/photons)]);
auto0 = auto;

clf;
tcspc = zeros(2^15, 2);
rate = [];
time = 0;
jj = 1;
cnt = 0;
tst = 1;

while tst>0
    [sync, micro, chan, bla, tst] = ht3v2read_old(name, [cnt+1, photons]);
    if tst>0
        cnt = cnt + tst;

        tcspc(:,1) = tcspc(:,1) + mHist(micro(chan==0),0:2^15-1);
        tcspc(:,2) = tcspc(:,2) + mHist(micro(chan==1),0:2^15-1);

        tim = sync*syncperiod;
        tim(chan==0) = tim(chan==0) + timeshift(micro(chan==0)+1,1);
        tim(chan==1) = tim(chan==1) + timeshift(micro(chan==1)+1,2);        
        time(jj) = tim(end) - tim(1);
        
        for k=1:ceil(automaxtime/150)
            tmp = tim(k+1:end)-tim(1:end-k);
            ind = chan(1:end-k)==0 & weight(micro(1:end-k)+1,1)==1 & chan(k+1:end)==1 & weight(micro(k+1:end)+1,2)==1;
            if sum(ind)>1
                auto(:,1,1,1,jj) = auto(:,1,1,1,jj) + mHist(tmp(ind),autotime);
            end
            ind = chan(1:end-k)==0 & weight(micro(1:end-k)+1,1)==1 & chan(k+1:end)==1 & weight(micro(k+1:end)+1,2)==2;
            if sum(ind)>1
                auto(:,1,2,1,jj) = auto(:,1,2,1,jj) + mHist(tmp(ind),autotime);
            end
            ind = chan(1:end-k)==0 & weight(micro(1:end-k)+1,1)==2 & chan(k+1:end)==1 & weight(micro(k+1:end)+1,2)==1;
            if sum(ind)>1
                auto(:,2,1,1,jj) = auto(:,2,1,1,jj) + mHist(tmp(ind),autotime);
            end
            ind = chan(1:end-k)==0 & weight(micro(1:end-k)+1,1)==2 & chan(k+1:end)==1 & weight(micro(k+1:end)+1,2)==2;
            if sum(ind)>1
                auto(:,2,2,1,jj) = auto(:,2,2,1,jj) + mHist(tmp(ind),autotime);
            end
            ind = chan(1:end-k)==1 & weight(micro(1:end-k)+1,2)==1 & chan(k+1:end)==0 & weight(micro(k+1:end)+1,1)==1;
            if sum(ind)>1
                auto(:,1,1,2,jj) = auto(:,1,1,2,jj) + mHist(tmp(ind),autotime);
            end
            ind = chan(1:end-k)==1 & weight(micro(1:end-k)+1,2)==1 & chan(k+1:end)==0 & weight(micro(k+1:end)+1,1)==2;
            if sum(ind)>1
                auto(:,1,2,2,jj) = auto(:,1,2,2,jj) + mHist(tmp(ind),autotime);
            end
            ind = chan(1:end-k)==1 & weight(micro(1:end-k)+1,2)==2 & chan(k+1:end)==0 & weight(micro(k+1:end)+1,1)==1;
            if sum(ind)>1
                auto(:,2,1,2,jj) = auto(:,2,1,2,jj) + mHist(tmp(ind),autotime);
            end
            ind = chan(1:end-k)==1 & weight(micro(1:end-k)+1,2)==2 & chan(k+1:end)==0 & weight(micro(k+1:end)+1,1)==2;
            if sum(ind)>1
                auto(:,2,2,2,jj) = auto(:,2,2,2,jj) + mHist(tmp(ind),autotime);
            end
        end

        rate(jj,1,1) = sum(chan==0 & weight(micro+1,1)==1)/time(jj);
        rate(jj,2,1) = sum(chan==0 & weight(micro+1,1)==2)/time(jj);
        rate(jj,1,2) = sum(chan==1 & weight(micro+1,1)==1)/time(jj);
        rate(jj,2,2) = sum(chan==1 & weight(micro+1,1)==2)/time(jj);

        subplot('position', [0.925 0.2 0.025 0.6]);
        bar(0,cnt/NCounts);
        axis([-0.4 0.4 0 1]);
        set(gca,'xtick',[],'ytick',[]);
        subplot('position', [0.1 0.1 0.7 0.8]);
        plot(-autotime(sum(auto(:,1,2,2,:),5)>0), sum(auto(sum(auto(:,1,2,2,:),5)>0,1,2,2,:),5), ...
            -autotime(sum(auto(:,2,1,2,:),5)>0), sum(auto(sum(auto(:,2,1,2,:),5)>0,2,1,2,:),5), ...
            autotime(sum(auto(:,1,2,1,:),5)>0), sum(auto(sum(auto(:,1,2,1,:),5)>0,1,2,1,:),5), ...
            autotime(sum(auto(:,2,1,1,:),5)>0), sum(auto(sum(auto(:,2,1,1,:),5)>0,2,1,1,:),5));
        drawnow;
    end
    jj = jj+1;
end
clf

tau = (0.5:32768)*microresolution;
tau(sum(tcspc,2)==0) = [];
tcspc(sum(tcspc,2)==0,:) = [];
% [ind,ind] = min(gradient(sum(tcspc,2)));
% tcspc = tcspc(1:ind-10,:);
% tau = tau(1:ind-10);
% [ind,ind] = max(tcspc(:,1));
% len = length(tau);
% tcspc(:,1) = tcspc(mod(len+(1:len)+ind-500,len)+1,1);
% tmp = [];
% for j=1:size(tcspc,1)
%     tmp(j) = tcspc(:,1)'*tcspc([j:end 1:j-1],2);
% end
% [tt, len] = mCluster(tmp>0.9*mean(tmp));
% ind = max(tmp(tt==1));
% num = 1:length(tmp);
% ind = num(tmp==ind);
% ind = ind(1);
% tcspc(:,2) = tcspc([ind:end 1:ind-1],2);
plot(tau,tcspc(:,1)/sum(tcspc(:,1)),tau,tcspc(:,2)/sum(tcspc(:,2)))

win = tau<2;
mm(1,:) = max(tcspc(win,:));
win = 4<tau & tau<6;
if not(isempty(max(tcspc(win,:))))
    mm(2,:) = max(tcspc(win,:));
end;

figure
plot(-autotime(sum(auto(:,1,2,2,:),5)>0), sum(auto(sum(auto(:,1,2,2,:),5)>0,1,2,2,:),5), ...
    -autotime(sum(auto(:,2,1,2,:),5)>0), sum(auto(sum(auto(:,2,1,2,:),5)>0,2,1,2,:),5), ...
    autotime(sum(auto(:,1,2,1,:),5)>0), sum(auto(sum(auto(:,1,2,1,:),5)>0,1,2,1,:),5), ...
    autotime(sum(auto(:,2,1,1,:),5)>0), sum(auto(sum(auto(:,2,1,1,:),5)>0,2,1,1,:),5));
xlabel('lag time (ns)'); ylabel('autocorrelation')

if 1
    t1 = autotime(sum(auto(:,1,2,1,:),5)>0);
    t2 = autotime(sum(auto(:,2,1,2,:),5)>0);
    a1 = squeeze(auto(sum(auto(:,1,2,1,:),5)>0,1,2,1,:));
    a2 = squeeze(auto(sum(auto(:,2,1,2,:),5)>0,2,1,2,:));
    p = Simplex('ExpFun',[100 1000],[],[],[],[], t1, [sum(a1,2), sum(a2,2)], 1);
    [err, c1, zz1, z1] = ExpFun(p, t1, sum(a1,2));
    z1 = ExpFun(p, c1, autotime);
    [err, c2, zz2, z2] = ExpFun(p, t2, sum(a2,2));
    z2 = ExpFun(p, c2, autotime);
    plot(-t2, sum(a2,2),'o', -autotime, z2, t1, sum(a1,2),'o', autotime,z1);
    ax = axis;
    axis([-automaxtime automaxtime ax(3) ax(4)])
    xlabel('time (ns)'); ylabel('autocorrelation (a.u.)')
    text(ax(1)+2*diff(ax(1:2))/3,ax(3)+2*diff(ax(3:4))/3,['\tau_{\itrot\rm} = ' mnum2str(p(1),3,1), ' ns']);
else
    p = Simplex('BSARotationFun',[0.1 50 2e3],[],[],[],[],autotime,auto);
    ax = axis;
    axis([-automaxtime automaxtime ax(3) ax(4)])
    xlabel('time (ns)'); ylabel('autocorrelation (a.u.)')
    text(ax(1)+0.51*diff(ax(1:2)),ax(3)+2*diff(ax(3:4))/3,['\tau_{\itrot\rm} = ' mnum2str(p(2),3,1), ' ns']);
end

tmp = [0 findstr(name,'\')];
name = [name(tmp(end)+1:end-4) 'All'];
eval(['save ''' name ''' bin tcspcraw tau tcspc mm autotime auto time rate microresolution syncperiod'])
eval(['print -dpng -r300 ' name])

return
beta = sum(auto,5);
t1a2 = flipud(autotime(beta(:,1,1,2)>0));
t1b2 = flipud(autotime(beta(:,2,2,2)>0));
a1a2 = flipud(beta(beta(:,1,1,2)>0,1,1,2));
a1b2 = flipud(beta(beta(:,2,2,2)>0,2,2,2));
t22 = flipud(autotime(beta(:,2,1,2)>0));
a22 = flipud(beta(beta(:,2,1,2)>0,2,1,2));
t32 = flipud(autotime(beta(:,1,2,2)>0));
a32 = flipud(beta(beta(:,1,2,2)>0,1,2,2));
t1a1 = autotime(beta(:,1,1,1)>0);
t1b1 = autotime(beta(:,2,2,1)>0);
a1a1 = beta(beta(:,1,1,1)>0,1,1,1);
a1b1 = beta(beta(:,2,2,1)>0,2,2,1);
t21 = autotime(beta(:,1,2,1)>0);
a21 = beta(beta(:,1,2,1)>0,1,2,1);
t31 = autotime(beta(:,2,1,1)>0);
a31 = beta(beta(:,2,1,1)>0,2,1,1);
plot(-t1a2,a1a2,t1a1,a1a1,-t1b2,a1b2,t1b1,a1b1,-t22,a22,t21,a21,-t32,a32,t31,a31)
xlabel('time (ns)'); ylabel('autocorrelation (a.u.)')

t = unique([t1a1; t1a2; t1b1; t1b2; t21; t22; t31; t32]);
t(t<min(t21)) = [];

a11i = interp1(t1a1,a1a1,t,'cubic');
a12i = interp1(t1a2,a1a2,t,'cubic');
a11i = a11i + interp1(t1b1,a1b1,t,'cubic');
a12i = a12i + interp1(t1b2,a1b2,t,'cubic');
a21i = interp1(t21,a21,t,'cubic');
a22i = interp1(t22,a22,t,'cubic');
a31i = interp1(t31,a31,t,'cubic');
a32i = interp1(t32,a32,t,'cubic');
plot(t,sqrt(a12i.*a11i)+sqrt(a22i.*a21i)+sqrt(a32i.*a31i))
p = Simplex('ExpFun',[20 1e3],[],[],[],[],t,sqrt(a12i.*a11i)+sqrt(a22i.*a21i)+sqrt(a32i.*a31i),1)



return % anisotropy
plot(cumsum(time),(sqrt(rate(:,1,1).*rate(:,2,2))-sqrt(rate(:,2,1).*rate(:,1,2)))./(sqrt(rate(:,1,1).*rate(:,2,2))+2*sqrt(rate(:,2,1).*rate(:,1,2))))
mean((sqrt(rate(:,1,1).*rate(:,2,2))-sqrt(rate(:,2,1).*rate(:,1,2)))./(sqrt(rate(:,1,1).*rate(:,2,2))+2*sqrt(rate(:,2,1).*rate(:,1,2))))

