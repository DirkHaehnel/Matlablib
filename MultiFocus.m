% function [res, head] = MultiFocus(name);

timewin = 0.1; % time window of data processing in s
photons = 1e6; % number of photons processed one at a time
smooth = 1; % smoothed autocorrelation
binwidth = 1e3; 
filt = 0; 
close all   

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

bin = bin';
tmpa = tcspcdata(:,1)-min(tcspcdata(tcspcdata(:,1)>0,1));
tmpb = tcspcdata(:,2)-min(tcspcdata(tcspcdata(:,2)>0,2));

win = (1:round(length(bin)/4))-30;

tcspc = []; tau = [];
for j=1:3
    [pos pos] = max(tmpa);
    ind = pos + win;
    tcspc = [tcspc tmpa(ind)];
    tau = [tau bin(ind)];
    tmpa(ind) = 0;
    [pos pos] = max(tmpb);
    ind = pos + win;
    tcspc = [tcspc tmpb(ind)];
    tau = [tau bin(ind)];
    tmpb(ind) = 0;
end
[ind,ind] = sort(tau(1,1:2:end));
ind = 2*(ind-1)+1;
tau(:,1:2:end) = tau(:,ind);
tau(:,2:2:end) = tau(:,ind+1);
tcspc(:,1:2:end) = tcspc(:,ind);
tcspc(:,2:2:end) = tcspc(:,ind+1);
semilogy(tau,tcspc); drawnow;

numrange = -20:20;
maxnum = length(numrange);
hnum = zeros(maxnum,maxnum,maxnum);
if filt
else
    time = 0;
    cnt = 0;
    tst = 1;
    while tst>0
        [tmpy, flv, tmpx, head, tst] = tttrRead(name,[cnt+1 photons]);        
        if tst>0
            cnt = cnt + tst;
            clear ft fb
            for p=1:3    
                ft(:,p) = flv==0 & tmpx>=tau(1,2*(p-1)+1) & tau(end,2*(p-1)+1)>=tmpx;
                fb(:,p) = flv==1 & tmpx>=tau(1,2*p) & tau(end,2*p)>=tmpx;                
            end
            dt = tmpy(end)-tmpy(1);
            time = time + dt*timeunit;
            nn1 = tttr2bin(tmpy, binwidth, ft(:,1) | ft(:,2) | ft(:,3)) - tttr2bin(tmpy, binwidth, fb(:,1) | fb(:,2) | fb(:,3));
            nn2 = tttr2bin(tmpy, binwidth, ft(:,1) | fb(:,1) | ft(:,2) | fb(:,2)) - tttr2bin(tmpy, binwidth, ft(:,3) | fb(:,3));
            nn3 = tttr2bin(tmpy, binwidth, ft(:,1) | fb(:,1) | ft(:,3) | fb(:,3)) - tttr2bin(tmpy, binwidth, ft(:,2) | fb(:,2));
            hnum = hnum + mHist3(nn1,nn2,nn3,-20:20,-20:20,-20:20);
                    
            subplot('position', [0.925 0.2 0.025 0.6]);
            bar(0,cnt/head.NCounts);
            axis([-0.4 0.4 0 1]);
            set(gca,'xtick',[],'ytick',[]);
            subplot('position', [0.1 0.1 0.7 0.8]);
            mim(log(sum(hnum,3))); 
            drawnow;
        end
    end
end

% res = struct('autotime', autotime);
% res = setfield(res, 'auto', auto/time);
% res = setfield(res, 'tau', head.Resolution*tau);
% res = setfield(res, 'tcspc', tcspc);

head.NCounts = cnt;

break

time = 0;
cnt = 0;
tst = 1;
kappa = [];
while tst>0
    [tmpy, flv, tmpx, head, tst] = tttrRead(name,[cnt+1 photons]);        
    if tst>0
        cnt = cnt + tst;
        clear ft fb
        for p=1:3    
            ft(:,p) = flv==0 & tmpx>=tau(1,2*(p-1)+1) & tau(end,p)>=tmpx;
            fb(:,p) = flv==1 & tmpx>=tau(1,2*p) & tau(end,p)>=tmpx;                
        end
        dt = tmpy(end)-tmpy(1);
        time = time + dt*timeunit;
        nn = tttr2bin(tmpy, binwidth, ft(:,1) | ft(:,2) | ft(:,3) | fb(:,1) | fb(:,2) | fb(:,3));
        nn1 = tttr2bin(tmpy, binwidth, ft(:,1) | ft(:,2) | ft(:,3)) - tttr2bin(tmpy, binwidth, fb(:,1) | fb(:,2) | fb(:,3));
        nn2 = tttr2bin(tmpy, binwidth, ft(:,1) | fb(:,1) | ft(:,2) | fb(:,2)) - tttr2bin(tmpy, binwidth, ft(:,3) | fb(:,3));
        nn3 = tttr2bin(tmpy, binwidth, ft(:,1) | fb(:,1) | ft(:,3) | fb(:,3)) - tttr2bin(tmpy, binwidth, ft(:,2) | fb(:,2));
        ind = nn>0;
        nn1 = nn1(ind);
        nn2 = nn2(ind);
        nn3 = nn3(ind);
        nn = nn(ind);        
        ind = nn1 + maxnum*(nn2-1) + maxnum^2*(nn3-1) - (1+maxnum+maxnum^2)*(numrange(1)-1);
        kappa = [kappa; nn./hnum(ind)];
        
        subplot('position', [0.925 0.2 0.025 0.6]);
        bar(0,cnt/head.NCounts);
        axis([-0.4 0.4 0 1]);
        set(gca,'xtick',[],'ytick',[]);
        subplot('position', [0.1 0.1 0.7 0.8]);
        mhist(kappa); 
        drawnow;
    end
end
