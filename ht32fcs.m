function [res, head] = ht32fcs(name, cnum, maxtime, cutoff, photons, timegate)

if nargin<6 || isempty(timegate)
    timegate = [];  % predefined timegates
end
if nargin<5 || isempty(photons)
    photons = 5e6;   % number of photons processed one at a time
end

close all

head = [];

if strcmp(name(end-2:end),'ht3')
	head = HT3_Read(name);
end

if ~isempty(head)
    Timeunit   = 1/(head.SyncRate);
    Resolution = head.Resolution*1e-9;
    NCounts    = head.Records;
else
    disp('Not a valid file!');
    return
end

if nargin<3 || isempty(maxtime)
    maxtime = 3;
end

if nargin<2 || isempty(cnum)
    cnum = 1;
end

Nsub = 10;

if length(maxtime)==2
    Nsub = maxtime(2);
end

Ncasc    = ceil(log2(maxtime(1)/Timeunit/Nsub+1));
autotime = Timeunit*cumsum(reshape(repmat(2.^[0 0:Ncasc-1],Nsub,1),(Ncasc+1)*Nsub,1));

tres = 15/Timeunit;

[bin, tcspcdata] = ReadTCSPC(name);

NChannels = numel(bin);

dnum = size(tcspcdata,2);
% dind = find(sum(tcspcdata,1)>0)-1;

%% Detect Laser-pulse Time-gates

if (cnum>1)&&(isempty(timegate))
    
    % Determine maxima of all detectors
    for k = 1:dnum        
        tmpdata(:,k) = smooth(tcspcdata(:,k),1e-10./Resolution);
    end
    
    irf_w = floor(0.9e-9/Resolution);
    win   = floor(NChannels/cnum);     % Max #channels / pulse;
  
    im = zeros(dnum, cnum);
    m  = zeros(dnum, cnum);

    [tpm,i1] = max(tmpdata);
    off      = floor(mod(i1,win)-irf_w);

    for l = 1:cnum
        for k = 1:dnum
          ind = off(k)+(l-1)*win+(1:win);
          ind(ind>NChannels) = [];
          [m(k,l),im(k,l)] = max(tmpdata(ind,k));
          im(k,l) = im(k,l) + off(k)+(l-1)*win;
        end
    end

    ind = m<(max(m,[],2)*ones(1,cnum)./100);
    
    for l = 1:cnum
        for k = 1:dnum
            if ind(k,l) == 1
                mi = find(im(k,:) == i1(k));
                im(k,l) = im(k,mi)+floor((l-mi)*NChannels/cnum);
            end
        end
    end
    
    %% Determinetime-gates for photon assignment

    Ngate = min([abs(im(:,1)-im(:,2)); NChannels-abs(im(:,1)-im(:,2))]);

    g = zeros(cnum*dnum, 4);
    
    for k = 1:dnum
        for l = 1:cnum
            d = im(k,l)-irf_w;
            if d > 0
                if d <= NChannels - Ngate
                    g(k+(l-1)*dnum,1:2) = [d d+Ngate-1];
                else
                    g(k+(l-1)*dnum,:) = [d NChannels 1 d+Ngate-NChannels-1];
                end;
            else
                d = d + NChannels;
                if d > NChannels - Ngate
                    g(k+(l-1)*dnum,:) = [d NChannels 1 d+Ngate-NChannels-1];
                else
                    g(k+(l-1)*dnum,1:2) = [d d+Ngate-1];
                end;
            end;
        end
    end

    channels = repmat(0:dnum-1,[1 cnum]);

    for n = 1:cnum*dnum
        mch = mean(g(n,1):(g(n,2)+g(n,3)));
        pie = floor(mch/(NChannels/(cnum-0.5)));
        timegate(1+channels(n)+pie*dnum,:) = g(n,:);
    end

elseif (cnum == 1)&&(isempty(timegate))
    timegate = repmat([1 NChannels 1 NChannels],[dnum 1]);
    channels = (0:dnum-1);
    Ngate    = NChannels;    

else 
    if size(timegate,2) == cnum
    elseif size(timegate,2) == cnum*dnum 
    end
    
    if timegate(2)>timegate(1)
        t1  = [timegate(1) timegate(2) 0 0];
    else
        t1  = [timegate(1) size(tcspcdata,1) 1 timegate(2)];
    end
    
    len = 1 + timegate(2) - timegate(1);
    if len<1
        len = len + Ngate;
    end
    
    if t1(3)== 0
        bin       = bin(t1(1):t1(2));
        tcspcdata = tcspcdata(t1(1):t1(2),:);
    else
        bin       = [bin(t1(1):t1(2)); bin(t1(3):t1(4))];
        tcspcdata = [tcspcdata(t1(1):t1(2), :); tcspcdata(t1(3):t1(4), :)];
    end
    
    t2 = 1+t1-t1(1);
    
end

%% Align Decay curves

tcspc = zeros(Ngate,cnum*dnum);
tau   = Resolution.*(1:Ngate);

for n = 1:cnum*dnum
    if (timegate(n,3)==0)
        tmp = tcspcdata(timegate(n,1):timegate(n,2),1+channels(n));
    else
        tmp = [tcspcdata(timegate(n,1):timegate(n,2),1+channels(n)); tcspcdata(timegate(n,3):timegate(n,4),1+channels(n))];
    end
    tcspc(:,n) = tmp;
end

if nargin>4 && ~isempty(cutoff)
    ind   = tau > cutoff;
    tau   = tau(ind);
    tcspc = tcspc(ind,:);
end

semilogy(Resolution*tau(:,:,1),tcspc(:,:,1)); drawnow

res.tau       = tau;
res.bin       = bin;
res.tcspcdata = tcspcdata;

auto = zeros(length(autotime),cnum*dnum,cnum*dnum);

rate  = [];
tcspc = [];
time = 0;
cnt  = 0;
cnt1 = 0;
num  = 1;
jj   = 1;
ovc  = 0;
y    = [];
x    = [];
ch   = [];
tmpy = [];
tmpx = [];
chan = [];

while num>0
    
    [ttmpy, ttmpx, tchan, markers, num, loc] = HT3_Read(name, [cnt+1 photons head.length]);
    
    cnt = cnt + num;
    
    if (num>0) && (~isempty(ttmpy))
        
        ttmpy(markers~=0) = [];
        ttmpx(markers~=0) = [];
        tchan(markers~=0) = [];

        if ~isempty(tmpy)
            ovc = ovc + tmpy(end);
        end
        
        tmpy = [tmpy; ttmpy+ovc];
        tmpx = [tmpx; ttmpx];
        chan = [chan; tchan];
        
        ovc = loc;
        
        while tmpy(end) > tres

            ind = find(tmpy>tres,1,'first');
            ind  = ind  - 1;
            cnt1 = cnt1 + ind;
            
            y  = tmpy(1:ind);
            x  = tmpx(1:ind);
            ch = chan(1:ind);

            tmpy(1:ind) = [];
            tmpx(1:ind) = [];
            chan(1:ind) = [];

            tmpy = tmpy -y(end);
            
            time(jj) = y(end)*Timeunit;

            flvp = zeros(length(y),dnum*cnum);

            for p = 1:dnum
                tmp = mHist(x(ch==p-1),bin);
                for q = 1:cnum
                    n = p+(q-1)*dnum;
                    if (timegate(n,3)==0)
                        tcspc(:,n,jj) = tmp(timegate(n,1):timegate(n,2));
                    else
                        tcspc(:,n,jj) = [tmp(timegate(n,1):timegate(n,2)); tmp(timegate(n,3):timegate(n,4))];
                    end                     
                end
            end
            
            for p=1:dnum*cnum
                flvp(:,p) = ch==channels(p) & ((x>=timegate(p,1) & timegate(p,2)>=x));
                if timegate(p,4) > 0
                    ind = ch==channels(p) & ((x>=timegate(p,3) & timegate(p,4)>=x));
                    flvp(ind,p) = 1;
                    y(ind) = y(ind)-1;
                end
            end

            rate(jj,:) = sum(abs(flvp))/time(jj);

            auto(:,:,:,jj) = tttr2xfcs(y, flvp, Ncasc, Nsub);

            subplot('position', [0.925 0.2 0.025 0.6]);
            bar(0,cnt1/NCounts);
            axis([-0.4 0.4 0 1]);
            set(gca,'xtick',[],'ytick',[]);
            subplot('position', [0.1 0.1 0.7 0.8]);
            semilogx(autotime, reshape(sum(auto,4),[size(auto,1) size(auto,2)*size(auto,3)])/sum(time));
            drawnow;
            jj = jj+1;
        end
    end
    
end

if ~isempty(tmpy)
    y  = tmpy(1:end);
    x  = tmpx(1:end);
    ch = chan(1:end);

    time(jj) = y(end)*Timeunit;

    flvp = zeros(length(y),dnum*cnum);

    for p = 1:dnum
        tmp = mHist(x(ch==p-1),bin);
        for q = 1:cnum
            n = p+(q-1)*dnum;
            if (timegate(n,3)==0)
                tcspc(:,n,jj) = tmp(timegate(n,1):timegate(n,2));
            else
                tcspc(:,n,jj) = [tmp(timegate(n,1):timegate(n,2)); tmp(timegate(n,3):timegate(n,4))];
            end
        end
    end

    for p=1:dnum*cnum
        flvp(:,p) = ch==channels(p) & ((x>=timegate(p,1) & timegate(p,2)>=x));
        if timegate(p,4) > 0
            ind = ch==channels(p) & ((x>=timegate(p,3) & timegate(p,4)>=x));
            flvp(ind,p) = 1;
            y(ind) = y(ind)-1;
        end
    end

    rate(jj,:) = sum(abs(flvp))/time(jj);

    auto(:,:,:,jj) = tttr2xfcs(y, flvp, Ncasc, Nsub);
end

for j=1:length(time)
    auto(:,:,:,j) = auto(:,:,:,j)/time(j)/Timeunit;
end

clf
semilogx(autotime, reshape(sum(auto,4),[size(auto,1) size(auto,2)*size(auto,3)])/sum(time));
res.tcspc    = tcspc;
res.autotime = autotime;
res.auto = auto;
res.rate = rate;
res.time = time;
head.NCounts = cnt;
