function [res, head] = MultiFocusFCS(name, order);

warning off MATLAB:DeprecatedLogicalAPI
set(0,'defaultFigurePosition',[100 100 660 520].*[0.2 0.5 1.5 1.2])
set(0,'defaultAxesFontName', 'Times', 'defaultAxesFontSize', 16, 'defaultLineLineWidth', 2, ...
    'defaultTextFontName', 'Times', 'defaultTextFontSize', 16);
set(0, 'defaultAxesColorOrder', [1 0 0; 0 0 1; 0 1 0; 0 1 1; 1 0 1; 1 1 0; 0.25 0.25 0.25]);

timewin = 0.1; % time window of data processing in s
photons = 1e6; % number of photons processed one at a time
smooth = 1; % smoothed autocorrelation
cut = 1; % cut autocorrelatin at 1 mus
filt = 0; 
close all   

[tmp, tmp, tmp, head] = tttrRead(name,[1 1]);
timeunit = head.GlobClock*1e-9;    

Nsub = 10;
Ncasc = 17;
autotime = cumsum(reshape(repmat(2.^[0 0:Ncasc-1],Nsub,1),(Ncasc+1)*Nsub,1));
autotime = autotime*timeunit;

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

win = (1:round(length(bin)/4))-30;

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

auto = zeros(length(autotime),6,6);
for p=1:1%6
    for q=1:1%6
        if filt
            % filter generation
            para = [tcspc(:,p) tcspc(:,q) ones(size(tcspc,1),2)];
            for j=1:size(para,2)
                para(:,j) = para(:,j)/sum(para(:,j));
            end % normalization
            mp = mod(p-1,2);
            mq = mod(q-1,2);
            weight = [tcspcdata(tau(1,p)<=bin & bin<=tau(end,p),mp+1) tcspcdata(tau(1,q)<=bin & bin<=tau(end,q),mq+1)];    
            master = (para(:,1:2:end)'*diag(1./weight(:,1))*para(:,1:2:end))\(diag(1./weight(:,1))*para(:,1:2:end))';
            master = [master; (para(:,2:2:end)'*diag(1./weight(:,2))*para(:,2:2:end))\(diag(1./weight(:,2))*para(:,2:2:end))']';
            
            time = 0;
            cnt = 0;
            tst = 1;
            while tst>0
                [tmpy, flv, tmpx, head, tst] = tttrRead(name,[cnt+1 photons]);        
                if tst>0
                    cnt = cnt + tst;
                    dt = tmpy(end)-tmpy(1);
                    time = time + dt*timeunit;
                    flvp = flv==mp & tmpx>=tau(1,p) & tau(end,p)>=tmpx;
                    flvq = flv==mq & tmpx>=tau(1,q) & tau(end,q)>=tmpx;                
                    ind = flvp | flvq;
                    tmpy = tmpy(ind);
                    flvp = flvp(ind);
                    flvq = flvq(ind);                
                    tmpx = tmpx(ind);            
                    tmpx(flvp) = tmpx(flvp) - tau(1,p) + 1;
                    if ~(p==q)
                        tmpx(flvq) = tmpx(flvq) - tau(1,q) + 1;            
                    end
                    flvp = double(flvp);
                    flvq = double(flvq);                
                    if nargin>1
                        auto(:,p,q) = auto(:,p,q) + dt*tttr2nfcs(tmpy, master(tmpx,1).*flvp, ...  
                            master(tmpx,3).*flvq, Ncasc, Nsub, order);
                    else
                        auto(:,p,q) = auto(:,p,q) + dt*tttr2fcs(tmpy, master(tmpx,1).*flvp, ...  
                            master(tmpx,3).*flvq, Ncasc, Nsub, smooth);
                    end                    
                    subplot('position', [0.925 0.2 0.025 0.6]);
                    bar(0,cnt/head.NCounts);
                    axis([-0.4 0.4 0 1]);
                    set(gca,'xtick',[],'ytick',[]);
                    subplot('position', [0.1 0.1 0.7 0.8]);
                    semilogx(autotime, auto(:,p,q)/time); 
                    drawnow;
                end
            end
        else
            mp = mod(p-1,2);
            mq = mod(q-1,2);
           
            time = 0;
            cnt = 0;
            tst = 1;
            while tst>0
                [tmpy, flv, tmpx, head, tst] = tttrRead(name,[cnt+1 photons]);        
                if tst>0
                    cnt = cnt + tst;
                    flvp = flv==mp & tmpx>=tau(1,p) & tau(end,p)>=tmpx;
                    flvq = flv==mq & tmpx>=tau(1,q) & tau(end,q)>=tmpx;                
                    dt = tmpy(end)-tmpy(1);
                    time = time + dt*timeunit;
                    flvp = double(flvp);
                    flvq = double(flvq);                
                    if nargin>1
                        auto(:,p,q) = auto(:,p,q) + dt*tttr2nfcs(tmpy, flvp, flvq, Ncasc, Nsub, 2);
                    else
                        auto(:,p,q) = auto(:,p,q) + dt*tttr2fcs(tmpy, flvp, flvq, Ncasc, Nsub, smooth);
                    end
                    subplot('position', [0.925 0.2 0.025 0.6]);
                    bar(0,cnt/head.NCounts);
                    axis([-0.4 0.4 0 1]);
                    set(gca,'xtick',[],'ytick',[]);
                    subplot('position', [0.1 0.1 0.7 0.8]);
                    semilogx(autotime, auto(:,p,q)/time); 
                    drawnow;
                end
            end
        end
    end
end
clf
auto = auto/time;
if cut
    ind = autotime<1e-6;
    auto(ind,:,:) = [];
    autotime(ind) = [];
end

set(0,'defaultAxesFontSize', 10);
for j=1:9
    subplot(3,3,j)
    switch j
        case 1
            p=1; q=2; s1='T1'; s2='B1';
        case 2
            p=3; q=4; s1='T2'; s2='B2';
        case 3
            p=5; q=6; s1='T3'; s2='B3';
        case 4
            p=1; q=3; s1='T1'; s2='T2';
        case 5
            p=3; q=5; s1='T2'; s2='T3';
        case 6
            p=5; q=1; s1='T3'; s2='T1';
        case 7
            p=2; q=4; s1='B1'; s2='B2';
        case 8
            p=4; q=6; s1='B2'; s2='B3';
        case 9
            p=6; q=2; s1='B3'; s2='B1';
    end        
    FCSScalePlot(autotime,[auto(:,p,p) auto(:,q,q) auto(:,p,q)+auto(:,q,p)],1)
    %semilogx(autotime,[auto(:,p,p)/auto(end,p,p) auto(:,q,q)/auto(end,q,q) (auto(:,p,q)+auto(:,q,p))/(auto(end,p,q)+auto(end,q,p))]-1)
    xlabel(''); ylabel(''); % axis([autotime(1) autotime(end) 0 1.5]);
    if mod(j,3)==1 ylabel('autocorrelation'); end
    if j>6 xlabel('time [s]'); end
    legend({[s1 '/' s1], [s2 '/' s2], [s1 '/' s2]});
end
set(0,'defaultAxesFontSize', 16);

res = setfield(res, 'autotime', autotime);
res = setfield(res, 'auto', auto/time);

head.NCounts = cnt;

warning on MATLAB:DeprecatedLogicalAPI

