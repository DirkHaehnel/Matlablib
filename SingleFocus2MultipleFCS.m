function [res, head] = SingleFocus2MultipleFCS(name, order, timebin);

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
Ncasc = ceil(log2(1/timeunit/timebin/Nsub+1)); % 1 sec max correlation time
autotime = cumsum(reshape(repmat(2.^[0 0:Ncasc-1],Nsub,1),(Ncasc+1)*Nsub,1));
autotime = autotime*timeunit*timebin;
auto = zeros(length(autotime),order,order,2,2,1);
intensity = zeros(order,2,1);
tcspc = zeros(head.NChannels,2,1);

time = [];
cnt = 0;
tst = 1;
rate = [];
indcnt = 1;
while tst>0
    [tmpy, flv, tmpx, head, tst] = tttrRead(name,[cnt+1 photons]);        
    if tst>0
        cnt = cnt + tst;
        dt = (tmpy(end)-tmpy(1))/timebin;
        time = [time dt];
        tcspc(:,1,indcnt) = mhist(tmpx(flv==0), 1:head.NChannels);
        tcspc(:,2,indcnt) = mhist(tmpx(flv==1), 1:head.NChannels);            
        for p=1:2
            flvp = flv==p-1;
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
        for p=1:2
            for q=(p+1):2
                mp = mod(p-1,2);
                mq = mod(q-1,2);
                flvp = flv==mp;
                flvq = flv==mq;                
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

res = struct('autotime', autotime);
res = setfield(res, 'auto', squeeze(auto));
res = setfield(res, 'intensity', intensity);
tmp = 1:head.NChannels;
res = setfield(res, 'tau', head.Resolution*tmp(sum(sum(tcspc,3),2)>0));
res = setfield(res, 'tcspc', tcspc(sum(sum(tcspc,3),2)>0,:,:));
res = setfield(res, 'time', time);
head.NCounts = cnt;

