function [res, head] = MultiFocus2FCS(name, cnum, maxtime, para, cutoff, photons, timegates)

if nargin<7 || isempty(timegates)
    timegates = [];  % predefined timegates
end
if nargin<6 || isempty(photons)
    photons = 1e6; % number of photons processed one at a time
end

close all

if strcmp(name(end-2:end),'pt3')
    head = pt3v2Read(name);
    Timeunit = head.Sync*1e-9;
    Resolution = head.Resolution;
    NChannels  = 2^12;
    NCounts   = head.Records;
elseif strcmp(name(end-2:end),'t3r')
    head = tttrRead(name);
    Timeunit = head.GlobClock*1e-9;
    Resolution = head.Resolution;
    NChannels  = 2^12;
    NCounts   = head.Records;
elseif strcmp(name(end-2:end),'ht3')
	head = ht3v2read(name);
    if isempty(head)
        if nargin<4 || isempty(para)
            Timeunit = 50e-9;
            Resolution = 2e-12;
        else
            Timeunit = para(1)*1e-9;
            Resolution = para(2)*1e-9;
        end
        NChannels  = ceil(Timeunit/Resolution);
        tmp = dir(name);
        NCounts = tmp.bytes/4;
    else
        Timeunit = 1/(head.SyncRate);
        Resolution = head.Resolution*1e-9;
        NChannels = 2^15;
        tmp = dir(name);
        NCounts = tmp.bytes/4;
    end
else
    disp('Not a valid file!');
    return
end

if nargin<3 || isempty(maxtime)
    maxtime = 3;
    Nsub = 10;
    Ncasc = ceil(log2(maxtime/Timeunit/Nsub)); % 3 sec max correlation time
else
    if length(maxtime)==2
        Nsub = maxtime(2);
    else
        Nsub = 10;
    end
    Ncasc = ceil(log2(maxtime(1)/Timeunit/Nsub));
end

%if 3e3*maxtime(1)*NCounts/head.AcquisitionTime<photons
%    photons = round(3e3*maxtime(1)*NCounts/head.AcquisitionTime);
%end
autotime = cumsum(reshape(repmat(2.^(0:Ncasc-1),Nsub,1),Ncasc*Nsub,1));
autotime = autotime*Timeunit;

if strcmp(name(end-2:end),'pt3')
    [tmpy, tmpx, flv, bla, tst] = pt3v2Read(name, [1, 1e5]);
elseif strcmp(name(end-2:end),'t3r')
    [tmpy, tmpx, flv, bla, tst] = tttrRead(name, [1, 1e5]);
elseif strcmp(name(end-2:end),'ht3')
    if isempty(head)
        [tmpy, tmpx, flv, special, tst, overcount, head] = ht3read(name, [1, 1e5]);
    else
        [tmpy, tmpx, flv, special, tst, head] = ht3v2read(name, [1, 1e5]);
    end
end

dind = unique(flv);
dnum = length(dind);
for j=1:length(dind)
    rev_dind(dind(j)+1) = j;
end
        
% Read TCSPC Data
[bin, tcspcdata] = ReadTCSPC(name);
tcspcdata2 = zeros(NChannels,dnum,ceil(head.Records/photons));

% Detect Laser Puls Time Gates
if isempty(timegates)
    [t1, len] = AutodetectTimeGates(tcspcdata, cnum);
else
    t1 = timegates(:,1);
    len = 1 + min(timegates(:,2) - timegates(:,1));
end

%t2 = 1+t1-t1(1);
%bin = [bin(t1:end); bin(1:t1(1)-1)];
%tcspcdata = [tcspcdata(t1:end,:); tcspcdata(1:t1(1)-1,:)];

for j=1:cnum
    for k=1:dnum
        tau(1:len,j,k) = bin(t1(j):t1(j)+len-1);
        tcspc(1:len,j,k) = tcspcdata(t1(j):t1(j)+len-1,k);
    end
end
if nargin>4 && ~isempty(cutoff)
    ind = tau(:,1,1) - tau(1,1,1) > cutoff/Resolution;
    tau = tau(ind,:,:);
    tcspc = tcspc(ind,:,:);
end
tcspc(any(any(tau==0,2),3),:,:)=[];
tau(any(any(tau==0,2),3),:,:)=[];
semilogy(Resolution*tau(:,:,1),tcspc(:,:,1)); drawnow
res.tau = Resolution*tau;
res.tcspc = tcspc;
res.bin = bin;
res.tcspcdata = tcspcdata;

auto = zeros(length(autotime),cnum*dnum,cnum*dnum);

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
    ind = flv>(cnum*dnum-1);
    tmpy(ind) = [];
    tmpx(ind) = [];
    flv(ind)  = [];
    
    cnt = cnt + tst;
    if (tst>0) && (~isempty(tmpy))
        dt = tmpy(end)-tmpy(1);
        time(jj) = dt*Timeunit;
        
        if (tst>0) && (~isempty(tmpy))
            for j=1:length(flv)
                dd = rev_dind(flv(j)+1);
                if (dd > 0)
                    tcspcdata2(tmpx(j)+1, dd, jj) = tcspcdata2(tmpx(j)+1, dd, jj) + 1;
                end
            end
        end        

        flvp = zeros(length(flv),cnum*dnum);
        r = 1;
        for p=1:dnum
            for q=1:cnum
                flvp(:,r) = flv==dind(p) & tmpx>=tau(1,q,p) & max(tau(:,q,p))>=tmpx;
                r = r+1;
            end
        end
        rate(jj,:) = sum(flvp)/dt/Timeunit;

        auto(:,:,:,jj) = tttr2xfcs(tmpy, flvp, Ncasc, Nsub);

        subplot('position', [0.925 0.2 0.025 0.6]);
        bar(0,cnt/NCounts);
        axis([-0.4 0.4 0 1]);
        set(gca,'xtick',[],'ytick',[]);
        subplot('position', [0.1 0.1 0.7 0.8]);
        semilogx(autotime, reshape(sum(auto,4),[size(auto,1) size(auto,2)*size(auto,3)])/sum(time));
        drawnow;
    end
    jj = jj+1;
end
for j=1:length(time)
    auto(:,:,:,j) = auto(:,:,:,j)/time(j)/Timeunit;
end

clf
semilogx(autotime, reshape(sum(auto,4),[size(auto,1) size(auto,2)*size(auto,3)])/sum(time));
res.autotime = autotime;

ind = sum(sum(tcspcdata2,3),2)==0;
tcspcdata2(ind,:,:) = [];
res.tcspcdata2 = tcspcdata2;

res.auto = auto;
res.rate = rate;
res.time = time;
head.NCounts = cnt;
