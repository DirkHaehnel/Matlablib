function [res, head] = FRET2FLCS(name, maxtime, mask, photons, taufcs)

if nargin<4 || isempty(photons)
    photons = 1e6; % number of photons processed one at a time
end

close all

if strcmp(name(end-2:end),'ht3')
	head = ht3v2read(name);
    Timeunit = 1/(head.SyncRate);
    Resolution = head.Resolution*1e-9;
    NChannels = 2^15;
    tmp = dir(name);
    NCounts = tmp.bytes/4;
    dind = [0 1 2 3]; % red red blue blue 
else
    disp('Not a valid file!');
    return
end

if nargin<2 || isempty(maxtime)
    Nsub = 10;
    Ncasc = ceil(log2(3/Timeunit/Nsub)); % 3 sec max correlation time
else
    if length(maxtime)==2
        Nsub = maxtime(2);
    else
        Nsub = 10;
    end
    Ncasc = ceil(log2(maxtime(1)/Timeunit/Nsub));
end
autotime = cumsum(reshape(repmat(2.^(0:Ncasc-1),Nsub,1),Ncasc*Nsub,1));
autotime = autotime*Timeunit;

tcspcdata = zeros(NChannels,2);
tst = 1;
cnt = 0;
bin = 0:NChannels-1;
h = waitbar(0,'Please wait...');
while tst>0
    if strcmp(name(end-2:end),'pt3')
        [tmpy, tmpx, flv, bla, tst] = pt3v2Read(name, [cnt+1, photons]);
    elseif strcmp(name(end-2:end),'t3r')
        [tmpy, tmpx, flv, bla, tst] = tttrRead(name, [cnt+1, photons]);
    elseif strcmp(name(end-2:end),'ht3')
        if isempty(head)
            [tmpy, tmpx, flv, special, tst, overcount, head] = ht3read(name, [cnt+1, photons]);
        else
            [tmpy, tmpx, flv, special, tst, head] = ht3v2read(name, [cnt+1, photons]);
        end
    end
    cnt = cnt+tst;
    waitbar(cnt/NCounts,h)
    if tst>0
        tcspcdata(:,1) = tcspcdata(:,1) + mHist(tmpx(flv==dind(1)),bin);
        tcspcdata(:,2) = tcspcdata(:,2) + mHist(tmpx(flv==dind(2)),bin);
    end
end
close(h)

ind = sum(tcspcdata,2)==0;
bin(ind) = [];
tcspcdata(ind,:) = [];
bin(1:30)=[];
tcspcdata(1:30,:) = [];
semilogy(bin,tcspcdata); drawnow;

bin = bin';
tmp = mean(tcspcdata,2);

tmp2 = sort(tmp);
baseline = mean(tmp2(1:ceil(size(tmp2)/5)));
peak = mean(tmp2(ceil(size(tmp2)*0.99):end));

[ind, num] = mCluster(tmp>(baseline+0.15*(peak-baseline)));
t1 = max([1, min(bin(ind==length(num)))-50]);
t2 = min(bin(ind==length(num)-1))-50;
if t1>t2
    tmp = t1; t1 = t2; t2 = tmp;
end
len = min([t2-t1,length(bin)-t2-20]);
tau(1:len,[1 3]) = bin(t1:t1+len-1)*[1 1];
tau(1:len,[2 4]) = bin(t2:t2+len-1)*[1 1];
tcspc(1:len,1) = tcspcdata(t1:t1+len-1, 1);
tcspc(1:len,3) = tcspcdata(t1:t1+len-1, 2);
tcspc(1:len,2) = tcspcdata(t2:t2+len-1, 1);
tcspc(1:len,4) = tcspcdata(t2:t2+len-1, 2);
if nargin>2 && ~isempty(cutoff)
    ind = tau(:,1) - tau(1,1) > cutoff/Resolution;
    tau = tau(ind,:);
    tcspc = tcspc(ind,:);
end
tcspc(any(tau==0,2),:)=[];
tau(any(tau==0,2),:)=[];
semilogy(Resolution*tau,tcspc); drawnow
res.tau = Resolution*tau;
res.tcspc = tcspc;
res.bin = bin;
res.tcspcdata = tcspcdata;

auto = zeros(length(autotime),4,4);
if nargin>4 && ~isempty(taufcs)
    autotau = auto;
end

rate = [];
time = 0;
cnt = 0;
tst = 1;
jj = 1;

while tst>0
    if strcmp(name(end-2:end),'pt3')
        [tmpy, tmpx, flv, bla, tst] = pt3v2Read(name, [cnt+1, photons]);
    elseif strcmp(name(end-2:end),'t3r')
        [tmpy, tmpx, flv, bla, tst] = tttrRead(name, [cnt+1, photons]);
    elseif strcmp(name(end-2:end),'ht3')
        if isempty(head)
            [tmpy, tmpx, flv, bla, tst] = ht3read(name, [cnt+1, photons]);
        else
            [tmpy, tmpx, flv, bla, tst] = ht3v2read(name, [cnt+1, photons]);
        end
    end
    if tst>0
        cnt = cnt + tst;
        dt = tmpy(end)-tmpy(1);
        time(jj) = dt*Timeunit;

        flvp = zeros(length(flv),4);
        for p=1:4
            flvp(:,p) = flv==dind(floor((p-1)/2)+1) & tmpx>=tau(1,p) & max(tau(:,p))>=tmpx;
        end
        rate(jj,:) = sum(flvp)/dt/Timeunit;

        auto(:,:,:,jj) = tttr2xfcs(tmpy, flvp, Ncasc, Nsub);
        if nargin>4 && ~isempty(taufcs)
            for p=1:4
                flvp(:,p) = flvp(:,p).*(tmpx-tau(1,p));
            end
            autotau(:,:,:,jj) = tttr2xfcs(tmpy, flvp, Ncasc, Nsub);
        end

        subplot('position', [0.925 0.2 0.025 0.6]);
        bar(0,cnt/NCounts);
        axis([-0.4 0.4 0 1]);
        set(gca,'xtick',[],'ytick',[]);
        subplot('position', [0.1 0.1 0.7 0.8]);
        semilogx(autotime, sum(auto(:,1,3,:)+auto(:,3,1,:)+auto(:,2,4,:)+auto(:,4,2,:),4)/sum(time),...
            autotime, sum(auto(:,1,4,:)+auto(:,4,1,:)+auto(:,2,3,:)+auto(:,3,2,:),4)/sum(time));
        drawnow;
    end
    jj = jj+1;
end
for j=1:length(time)
    auto(:,:,:,j) = auto(:,:,:,j)/time(j)/Timeunit;
    if nargin>5 && ~isempty(taufcs)
        autotau(:,:,:,j) = autotau(:,:,:,j)/time(j)/Timeunit;
    end
end

clf
if ~strcmp(name(end-2:end),'ht3')
    ind = autotime<5e-7;
    auto(ind,:,:,:) = [];
    autotime(ind) = [];
    if nargin>5 && ~isempty(taufcs)
        autotau(ind,:,:,:) = [];
    end
end

clf
semilogx(autotime, sum(auto(:,1,3,:)+auto(:,3,1,:)+auto(:,2,4,:)+auto(:,4,2,:),4),...
    autotime, sum(auto(:,1,4,:)+auto(:,4,1,:)+auto(:,2,3,:)+auto(:,3,2,:),4));
res.autotime = autotime;
res.auto = auto;
if nargin>5 && ~isempty(taufcs)
    res.autotau = autotau;
end
res.rate = rate;
res.time = time;
head.NCounts = cnt;

return

% calculated tau correlation function
ind = res.autotime<1e-2 & res.autotime>5e-6;
clear z
for j = 1:4 
    z(:,j) = abs((mean(res.autotau(ind,j,j,:),4)-mean(mean(res.autotau(end-10:end,j,j,:),4)))./(mean(res.auto(ind,j,j,:),4)-mean(mean(res.auto(end-10:end,j,j,:),4))))*head.Resolution^2;
end

% plot infinite limit estimate
semilogx(res.autotime,mean(res.auto(:,j,j,:),4),res.autotime,mean(mean(res.auto(end-10:end,j,j,:),4))*ones(size(res.autotime)))

% fit tau correlation
j=3; p=Simplex('ExpFun',[1e-5 1e-3],[],[],[],[],res.autotime(ind),z(:,j),3)

% test
ind = res.autotime<1e-1 & res.autotime>5e-6;
clear z
[y,t] = FCSCrossRead(res);
tmp = res; tmp.auto=tmp.autotau;
[yy,t] = FCSCrossRead(tmp);
for j=1:4
    z(:,j) = abs((mean(yy(ind,j,:),3)-mean(mean(yy(end-10:end,j,:),3)))./(mean(y(ind,j,:),3)-mean(mean(y(end-10:end,j,:),3))))*head.Resolution^2;
end
