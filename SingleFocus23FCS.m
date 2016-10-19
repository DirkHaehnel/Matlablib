function [res, head] = SingleFocus23FCS(name, maxtime, para, cutoff, photons)

if nargin<5 || isempty(photons)
    photons = 1e6; % number of photons processed one at a time
end

close all

if strcmp(name(end-2:end),'pt3')
    head = pt3v2Read(name);
    Timeunit = head.Sync*1e-9;
    Resolution = head.Resolution;
    NChannels  = 2^12;
    NCounts   = head.Records;
    dind = [1 3];
elseif strcmp(name(end-2:end),'t3r')
    head = tttrRead(name);
    Timeunit = head.GlobClock*1e-9;
    Resolution = head.Resolution;
    NChannels  = 2^12;
    NCounts   = head.Records;
	dind = [0 1];
elseif strcmp(name(end-2:end),'ht3')
    if nargin<3 || isempty(para)
        Timeunit = 50e-9;
        Resolution = 2e-12;
    else
        Timeunit = para(1)*1e-9;
        Resolution = para(2)*1e-9;
    end
    NChannels  = ceil(Timeunit/Resolution);
    tmp = dir(name);
    NCounts = tmp.bytes/4;
    dind = [0 1];
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
while tst>0
    if strcmp(name(end-2:end),'pt3')
        [tmpy, tmpx, flv, bla, tst] = pt3v2Read(name, [cnt+1, photons]);
    elseif strcmp(name(end-2:end),'t3r')
        [tmpy, tmpx, flv, bla, tst] = tttrRead(name, [cnt+1, photons]);
    elseif strcmp(name(end-2:end),'ht3')
        [tmpy, tmpx, flv, special, tst, overcount, head] = ht3read(name, [cnt+1, photons]);
    end
    cnt = cnt+tst
    if tst>0
        tcspcdata(:,1) = tcspcdata(:,1) + mHist(tmpx(flv==dind(1)),bin);
        tcspcdata(:,2) = tcspcdata(:,2) + mHist(tmpx(flv==dind(2)),bin);
    end
end

if nargin>3 && ~isempty(cutoff)
    ind = bin - bin(1) > cutoff/Resolution;
    bin = bin(ind);
    tcspcdata = tcspcdata(ind,:);
end

ind = sum(tcspcdata,2)==0;
bin(ind) = [];
tcspcdata(ind,:) = [];
bin(1:30)=[];
tcspcdata(1:30,:) = [];
semilogy(bin,tcspcdata); drawnow;

res.bin = bin;
res.tcspcdata = tcspcdata;

auto3 = zeros(length(autotime),2,2,2);
auto2 = zeros(length(autotime),2,2);
auto = auto2;

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
        [tmpy, tmpx, flv, bla, tst] = ht3read(name, [cnt+1, photons]);
    end
    if tst>0
        cnt = cnt + tst;
        dt = tmpy(end)-tmpy(1);
        time(jj) = dt*Timeunit;

        flvp = zeros(length(flv),2);
        for p=1:2
            flvp(:,p) = flv==dind(p);
        end
        rate(jj,:) = sum(flvp)/dt/Timeunit;

        [auto3(:,:,:,:,jj) auto2(:,:,:,jj) auto(:,:,:,jj)] = tttr23fcs(tmpy, flvp, Ncasc, Nsub);

        subplot('position', [0.925 0.2 0.025 0.6]);
        bar(0,cnt/NCounts);
        axis([-0.4 0.4 0 1]);
        set(gca,'xtick',[],'ytick',[]);
        subplot('position', [0.1 0.1 0.7 0.8]);
        semilogx(autotime,sum(auto3(:,1,1,1,:),5))
        drawnow;
    end
    jj = jj+1;
end
for j=1:length(time)
    auto3(:,:,:,:,j) = auto3(:,:,:,:,j)/time(j)/Timeunit^2;
    auto2(:,:,:,j) = auto2(:,:,:,j)/time(j)/Timeunit;
    auto(:,:,:,j) = auto(:,:,:,j)/time(j)/Timeunit;
end

clf
if ~strcmp(name(end-2:end),'ht3')
    ind = autotime<5e-7;
    auto3(ind,:,:,:,:) = [];
    auto2(ind,:,:,:) = [];
    auto(ind,:,:,:) = [];
    autotime(ind) = [];
end

clf
semilogx(autotime,sum(auto3(:,1,1,1,:),5))
res.autotime = autotime;
res.auto3 = auto3;
res.auto2 = auto2;
res.auto = auto;
res.rate = rate;
res.time = time;
head.NCounts = cnt;

