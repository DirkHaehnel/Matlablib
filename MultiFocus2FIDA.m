function [res, head] = MultiFocus2FIDA(name, timebin);

photons = 1e6; % number of photons processed one at a time

close all   

if nargin<2 | isempty(timebin)
    timebin = 1e3; % 100 microsseconds
end

[tmp, tmp, tmp, head] = tttrRead(name,[1 1]);
timeunit = head.GlobClock*1e-9;    

tcspcdata = zeros(head.NChannels,2);
tst = 1;
cnt = 0;
bin = 1:head.NChannels;
while tst>0
    [tmpy, flv, tmpx, head, tst] = tttrRead(name,[cnt+1 photons]);  
    cnt = cnt+tst;
    if tst>0
        tcspcdata(:,1) = tcspcdata(:,1) + mhist(tmpx(flv==0),bin);
        tcspcdata(:,2) = tcspcdata(:,2) + mhist(tmpx(flv==1),bin);
    end
end
ind = sum(tcspcdata,2)==0;
bin(ind) = [];
tcspcdata(ind,:) = [];
bin(1:30)=[];
tcspcdata(1:30,:) = [];
semilogy(bin,tcspcdata); drawnow;

tmpa = tcspcdata(:,1); %-min(tcspcdata(tcspcdata(:,1)>0,1));
tmpb = tcspcdata(:,2); %-min(tcspcdata(tcspcdata(:,2)>0,2));
bin = bin';

win = (1:round(length(bin)/3.5))-30;

tcspc = zeros(length(win),6); tau = tcspc;
for j=1:3
    [pos pos] = max(tmpa);
    ind = pos + win;
    if ind(1)<=0 
        win(1:sum(ind<=0)) = [];
        tcspc(1:sum(ind<=0),:) = [];
        tau(1:sum(ind<=0),:) = [];        
        ind = pos + win;
    end
    if ind(end)>length(tmpa) 
        win(end-sum(ind>length(tmpa)):end) = [];
        tcspc(end-sum(ind>length(tmpa)):end,:) = [];
        tau(end-sum(ind>length(tmpa)):end,:) = [];
        ind = pos + win;
    end
    tcspc(:,2*j-1) = tmpa(ind);
    tau(:,2*j-1) = bin(ind);
    tmpa(ind) = 0;
    [pos pos] = max(tmpb);
    ind = pos + win;
    if ind(1)<=0 
        win(1:sum(ind<=0)) = [];
        tcspc(1:sum(ind<=0),:) = [];
        tau(1:sum(ind<=0),:) = [];
        ind = pos + win;
    end
    if ind(end)>length(tmpb) 
        win(end-sum(ind>length(tmpb)):end) = [];
        tcspc(end-sum(ind>length(tmpb)):end,:) = [];
        tau(end-sum(ind>length(tmpb)):end,:) = [];
        ind = pos + win;
    end
    tcspc(:,2*j) = tmpb(ind);
    tau(:,2*j) = bin(ind);
    tmpb(ind) = 0;
end
[ind,ind] = sort(tau(1,1:2:end));
ind = 2*(ind-1)+1;
tau(:,1:2:end) = tau(:,ind);
tcspc(:,1:2:end) = tcspc(:,ind);
[ind,ind] = sort(tau(1,2:2:end));
ind = 2*ind;
tau(:,2:2:end) = tau(:,ind);
tcspc(:,2:2:end) = tcspc(:,ind);
tcspc(:,1:2:end) = tcspc(:,1:2:end)-min(min(tcspc(:,1:2:end)));
tcspc(:,2:2:end) = tcspc(:,2:2:end)-min(min(tcspc(:,2:2:end)));
semilogy(tau,tcspc); drawnow
res = struct('tau', head.Resolution*tau);
res = setfield(res, 'tcspc', tcspc);

fida = [];

time = 0;
cnt = 0;
tst = 1;
while tst>0
    [tmpy, flv, tmpx, head, tst] = tttrRead(name,[cnt+1 photons]);        
    if tst>0
        cnt = cnt + tst;
        dt = tmpy(end)-tmpy(1);
        time = time + dt*timeunit;
        flag = zeros(size(flv));
        for p=1:6
            mp = mod(p-1,2);
            flag(flv==mp & tmpx>=tau(1,p) & tau(end,p)>=tmpx) = p;
        end
        tmp = tttr2bin(tmpy, timebin, flag);
        tmp = tmp(:,2:end);
        tmp(prod(tmp,2)==0,:) = [];
        fida = [fida tmp'];
            
        subplot('position', [0.925 0.2 0.025 0.6]);
        bar(0,cnt/head.NCounts);
        axis([-0.4 0.4 0 1]);
        set(gca,'xtick',[],'ytick',[]);
        subplot('position', [0.1 0.1 0.7 0.8]);
        mhist(sum(fida),1:max(sum(fida))); 
        drawnow;
    end
end
clf
mhist(sum(fida),1:max(sum(fida))); 

res = setfield(res, 'fida', fida);
head.NCounts = cnt;

