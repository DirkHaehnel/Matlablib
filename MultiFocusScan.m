function [tag, tim, head, tau, tcspc] = MultiFocusScan(name);

close all
set(0,'defaultFigurePosition',[100 100 660 520].*[0.2 0.5 1.5 1.2])
set(0,'defaultAxesFontName', 'Times', 'defaultAxesFontSize', 16, 'defaultLineLineWidth', 2, ...
    'defaultTextFontName', 'Times', 'defaultTextFontSize', 16);
set(0, 'defaultAxesColorOrder', [1 0 0; 0 0 1; 0 1 0; 0 1 1; 1 0 1; 1 1 0; 0.25 0.25 0.25]);

photons = 1e6; % number of photons processed one at a time

[tmp, tmp, tmp, head] = tttrRead(name,[0 0]);
data = zeros(head.NChannels,2);
tst = 1;
cnt = 0;
bin = 1:head.NChannels;
while tst>0 & cnt<1e6
    [tmpy, flv, tmpx, head, tst] = tttrRead(name,[cnt+1 photons]);  
    cnt = cnt+tst;
    if tst>0
        data(:,1) = data(:,1) + mhist(tmpx(flv==0),bin);
        data(:,2) = data(:,2) + mhist(tmpx(flv==1),bin);
    end
end
ind = sum(data,2)==0;
bin(ind) = [];
tcspcdata = data;
tcspcdata(ind,:) = [];
bin(1:30)=[];
tcspcdata(1:30,:) = [];

tmpa = tcspcdata(:,1)-min(tcspcdata(tcspcdata(:,1)>0,1));
tmpb = tcspcdata(:,2)-min(tcspcdata(tcspcdata(:,2)>0,2));
bin = bin';

win = (1:round(length(bin)/4))-30;

tcspc = zeros(length(win),6); tau = tcspc;
for j=1:3
    [pos pos] = max(tmpa);
    ind = pos + win;
    if ind(1)<=0 
        win(1:sum(ind<=0)) = [];
        tcspc(1:sum(ind<=0),:) = 0;
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

filter = zeros(head.NChannels,6);
col = ones(size(tau,1),1);
for j=1:3
    filter(tau(:,2*j-1),2*j-1) = col;
    filter(tau(:,2*j),2*j) = col;    
end

[tag, tim] = PIE710Read(name, filter);

for j=1:3 
    for k=1:2 
        subplot(2,3,j+3*(k-1)); 
        mim(tim(:,:,j,k)); 
        if j==1
            axis on
            set(gca,'xticklabel',[],'yticklabel',[]);
            ylabel(['SPAD # ' int2str(k)]);
        end
        if k==1
            title(['Laser # ' int2str(j)]);
        end
    end
end

