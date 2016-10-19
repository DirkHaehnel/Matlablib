function [y, x, flag, tcspc, tau, head] = MultiFocusRead(name, num);

[tmp, tmp, tmp, head] = tttrRead(name,[1 1]);
bin = 1:head.NChannels;
tcspcdata = zeros(head.NChannels,2);
[y, flv, x, head, tst] = tttrRead(name,num);  
tcspcdata(:,1) = tcspcdata(:,1) + mhist(x(flv==0),bin);
tcspcdata(:,2) = tcspcdata(:,2) + mhist(x(flv==1),bin);
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

flag = zeros(size(y));
for p=1:6
    mp = mod(p-1,2);
    flag(flv==mp & x>=tau(1,p) & tau(end,p)>=x) = p;
end
y(flag==0) = [];
x(flag==0) = [];
flag(flag==0) = [];
flag = flag-1;
