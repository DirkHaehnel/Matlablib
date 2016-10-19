function [res, head] = OneFocus2BlindFCS(name, maxtime, autobin, filt)

photons = 1e6; % number of photons processed one at a time

close all

if strcmp(name(end-2:end),'pt3')
    head = pt3v2read(name);
    Timeunit = head.Sync*1e-9;
    Resolution = head.Resolution*1e-9;
    NCounts   = head.Records;
    dind = [1 3];
elseif strcmp(name(end-2:end),'t3r')
    head = tttrRead(name);
    Timeunit = head.GlobClock*1e-9;
    Resolution = head.Resolution*1e-9;
    NCounts   = head.Records;
	dind = [0 1];
else
    disp('Not a valid file!');
    return
end
NChannels  = 2^12;

if nargin<3 || isempty(autobin)
    autobin = 1;
else
    Timeunit = autobin*Timeunit;
end
if nargin<2 || isempty(maxtime)
    Nsub = 10;
    Ncasc = ceil(log2(3/Resolution/Nsub)); % 3 sec max correlation time
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

tcspc= zeros(NChannels,2);
tst = 1;
cnt = 0;
bin = (1:NChannels)';
while tst>0
    if strcmp(name(end-2:end),'pt3')
        [tmpy, tmpx, flv, bla, tst] = pt3v2read(name, [cnt+1, photons]);
    else
        [tmpy, tmpx, flv, bla, tst] = tttrread(name, [cnt+1, photons]);
    end
    cnt = cnt+tst;
    if tst>0
        tcspc(:,1) = tcspc(:,1) + mHist(tmpx(flv==dind(1)),bin);
        tcspc(:,2) = tcspc(:,2) + mHist(tmpx(flv==dind(2)),bin);
    end
    if nargin<4 || isempty(filt)
        tst = 0;
    end
end
ind = sum(tcspc,2)==0;
bin(ind) = [];
tcspc(ind,:) = [];
bin(1:30)=[];
tcspc(1:30,:) = [];
semilogy(bin,tcspc); drawnow;

res.bin = bin;
res.tcspc = tcspc;

if nargin>3 && ~isempty(filt)
    tau1 = filt(:,1);
    tau2 = filt(:,2);
    filt = filt(:,3:end);
    len = size(filt,2)/2;
    
    % filter generation
    % para = [dye1 channel#1; dye1 channel#2; dye2 channel#1; dye2 channel#2];
    for j=1:2*len
        filt(:,j) = filt(:,j)/sum(filt(:,j));
    end % normalization
    weight = [tcspc(tau1(1)<=bin & bin<=tau1(end),1) tcspc(tau2(1)<=bin & bin<=tau2(end),2)];    
    master = (filt(:,1:2:end)'*diag(1./weight(:,1))*filt(:,1:2:end))\(diag(1./weight(:,1))*filt(:,1:2:end))';
    master = [master; (filt(:,2:2:end)'*diag(1./weight(:,2))*filt(:,2:2:end))\(diag(1./weight(:,2))*filt(:,2:2:end))']';
    
    auto = zeros(length(autotime),2*len,2*len);

    time = 0;
    cnt = 0;
    tst = 1;
    jj = 1;
    while tst>0
        if strcmp(name(end-2:end),'pt3')
            [tmpy, tmpx, flv, bla, tst] = pt3v2read(name, [cnt+1, photons]);
        else
            [tmpy, tmpx, flv, bla, tst] = tttrread(name, [cnt+1, photons]);
        end
        if tst>0
            cnt = cnt + tst;
            dt = tmpy(end)-tmpy(1);
            time(jj) = dt*Timeunit;

            ind = (flv==1 & (tmpx<tau1(1) | tmpx>tau1(end))) | (flv==3 & (tmpx<tau2(1) | tmpx>tau2(end)));
            tmpy(ind) = [];
            flv(ind) = [];
            tmpx(ind) = [];
            
            flvp = zeros(length(flv),size(master,2));
            for p=1:len
                flvp(flv==1,p) = master(tmpx(flv==1)-tau1(1)+1,p);
                flvp(flv==3,len+p) = master(tmpx(flv==3)-tau2(1)+1,len+p);
            end
            rate(jj,:) = sum(flvp)/dt/Timeunit;
            
            auto(:,:,:,jj) = tttr2xfcs(round(tmpy/autobin), flvp, Ncasc, Nsub);

            subplot('position', [0.92 0.2 0.025 0.6]);
            bar(0,cnt/NCounts);
            axis([-0.4 0.4 0 1]);
            set(gca,'xtick',[],'ytick',[]);
            subplot('position', [0.1 0.1 0.7 0.8]);
            semilogx(autotime(1:end), squeeze(sum(auto(1:end,1,3,:)+auto(1:end,3,1,:),4))); 
            drawnow;
        end
        jj = jj + 1;
    end
    for j=1:length(time)
        auto(:,:,:,j) = auto(:,:,:,j)/time(j)/Timeunit;
    end
    res.master = master;
else
    auto = zeros(length(autotime),2,2);

    rate = [];
    time = 0;
    cnt = 0;
    tst = 1;
    jj = 1;
    while tst>0
        if strcmp(name(end-2:end),'pt3')
            [tmpy, tmpx, flv, bla, tst] = pt3v2read(name, [cnt+1, photons]);
        else
            [tmpy, tmpx, flv, bla, tst] = tttrread(name, [cnt+1, photons]);
        end
        if tst>0
            cnt = cnt + tst;
            dt = tmpy(end)-tmpy(1);
            time(jj) = dt*Timeunit;

            flvp = zeros(length(flv),2);
            for p=1:2
                flvp(:,p) = flv==(2*p-1);
            end
            rate(jj,:) = sum(flvp)/dt/Timeunit;

            auto(:,:,:,jj) = tttr2xfcs(round(tmpy/autobin), flvp, Ncasc, Nsub);

            subplot('position', [0.925 0.2 0.025 0.6]);
            bar(0,cnt/NCounts);
            axis([-0.4 0.4 0 1]);
            set(gca,'xtick',[],'ytick',[]);
            subplot('position', [0.1 0.1 0.7 0.8]);
            semilogx(autotime, sum(auto(:,1,2,:)+auto(:,2,1,:),4));
            drawnow;
        end
        jj = jj+1;
    end
    for j=1:length(time)
        auto(:,:,:,j) = auto(:,:,:,j)/time(j)/Timeunit;
    end
end

clf
ind = autotime<5e-7;
auto(ind,:,:,:) = [];
autotime(ind) = [];

clf
semilogx(autotime, sum(auto(:,1,2,:)+auto(:,2,1,:),4));
res.autotime = autotime;
res.auto = auto;
res.rate = rate;
res.time = time;
head.NCounts = cnt;


