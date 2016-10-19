function [res, head] = MultiFocus2MultipleFCS(name, order, timebin);

photons = 1e6; % number of photons processed one at a time

close all   

if nargin<2 | isempty(order)
    order=4; 
end
if nargin<3 | isempty(timebin)
    timebin=100; 
end

[tmp, tmp, tmp, head] = tttrRead(name,[1 1]);
timeunit = head.GlobClock*1e-9;    

Nsub = 10;
Ncasc = ceil(log2(5/timeunit/timebin/Nsub+1)); % 5 sec max correlation time
autotime = cumsum(reshape(repmat(2.^[0 0:Ncasc-1],Nsub,1),(Ncasc+1)*Nsub,1));
autotime = autotime*timeunit*timebin;

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

auto = zeros(length(autotime),order,order,6,6,1);
intensity = zeros(order,6,1);

time = [];
cnt = 0;
tst = 1;
indcnt = 1;
while tst>0
    [tmpy, flv, tmpx, head, tst] = tttrRead(name,[cnt+1 photons]);        
    if tst>0
        cnt = cnt + tst;
        dt = (tmpy(end)-tmpy(1))/timebin;
        time = [time dt];
        for p=1:6
            mp = mod(p-1,2);
            flvp = flv==mp & tmpx>=tau(1,p) & tau(end,p)>=tmpx;
            tmp = ceil(tmpy(flvp)/timebin); flvp = double(flvp(flvp)); 
            t = 1:length(tmp);
            dtmp = diff(tmp);
            t1 = t([1; dtmp]>0 & [dtmp; 1]==0);
            t2 = t([1; dtmp]==0 & [dtmp; 1]>0);
            for j=1:length(t1)
                flvp(t1(j)) = sum(flvp(t1(j):t2(j)));
            end
            ind = [1; dtmp]==0;
            tmp(ind) = [];
            flvp(ind) = [];
            
            for r=1:order
                intensity(r,p,indcnt) = sum(flvp.^r)/dt;
                for s=1:order
                    auto(:,r,s,p,p,indcnt) = tttr2fcs(tmp, flvp.^r, flvp.^s, Ncasc, Nsub, 1)/dt;
                end
            end
        end
        for p=1:6
            for q=(p+1):6
                mp = mod(p-1,2);
                mq = mod(q-1,2);
                flvp = flv==mp & tmpx>=tau(1,p) & tau(end,p)>=tmpx;
                flvq = flv==mq & tmpx>=tau(1,q) & tau(end,q)>=tmpx;                
                ind = flvp | flvq;
                tmp = ceil(tmpy(ind)/timebin); flvp = double(flvp(ind)); flvq = double(flvq(ind));
                t = 1:length(tmp);
                dtmp = diff(tmp);
                t1 = t([1; dtmp]>0 & [dtmp; 1]==0);
                t2 = t([1; dtmp]==0 & [dtmp; 1]>0);
                for j=1:length(t1)
                    flvp(t1(j)) = sum(flvp(t1(j):t2(j)));
                    flvq(t1(j)) = sum(flvq(t1(j):t2(j)));
                end
                ind = [1; dtmp]==0;
                tmp(ind) = [];
                flvp(ind) = [];
                flvq(ind) = [];
                
                for r=1:order
                    for s=1:order
                        auto(:,r,s,p,q,indcnt) = tttr2fcs(tmp, flvp.^r, flvq.^s, Ncasc, Nsub, 1)/dt;
                        auto(:,r,s,q,p,indcnt) = tttr2fcs(tmp, flvq.^r, flvp.^s, Ncasc, Nsub, 1)/dt;
                    end
                end

                subplot('position', [0.925 0.2 0.025 0.6]);
                bar(0,cnt/head.NCounts);
                axis([-0.4 0.4 0 1]);
                set(gca,'xtick',[],'ytick',[]);
                subplot('position', [0.1 0.1 0.7 0.8]);
                semilogx(autotime, sum(auto(:,1,1,p,q,:),6)); 
                drawnow;
            end
        end
    end
    indcnt = indcnt + 1;
end
clf

res = setfield(res, 'autotime', autotime);
res = setfield(res, 'auto', squeeze(auto));
res = setfield(res, 'intensity', intensity);
res = setfield(res, 'time', time);
head.NCounts = cnt;

