function [res, head] = SingleFocusCW2FCS_GridScan(name, maxtime, photons)

% Parameters
% name :        Filename

head = ht3v2read_head(name);

photons = 1e6; % number of photons processed one at a time
reduction = 10;    % online binning factor

isMTSGrid = (regexp(head.CommentField, 'MT Scan Grid') == 1);

if (isMTSGrid)
    Timeunit = head.Resolution;                 % timeunit in ns
    Timeunit2 = 1/(head.SyncRate);              % timeunit in s
    microresolution = head.Resolution;          % resolution in ns
    syncperiod = (1 / head.SyncRate) * 1e9;     % syncperiod in ns

    if nargin<2 || isempty(maxtime)
        maxtime = 1e3;
        Nsub = 10;
        Ncasc = ceil(log2(1e3/Timeunit/Nsub));    % 1 µsec max correlation time
        Ncasc2 = ceil(log2(0.1/Timeunit2/Nsub));    % 0.1 sec max correlation time
    else
        if length(maxtime)==2
            Nsub = maxtime(2);
        else
            Nsub = 10;
        end
        Ncasc = ceil(log2(maxtime(1)/Timeunit/Nsub));
        Ncasc2 = ceil(log2(0.1/Timeunit2/Nsub));    % 0.1 sec max correlation time
    end
    autotime = cumsum(reshape(repmat(2.^(0:Ncasc-1),Nsub,1),Ncasc*Nsub,1));
    autotime = autotime*Timeunit;      % in ns

    autotime2 = cumsum(reshape(repmat(2.^(0:Ncasc2-1),Nsub,1),Ncasc2*Nsub,1));
    autotime2 = autotime2*Timeunit2;   % in s

    autowidth = sum(horzcat(vertcat(0, diff(autotime)/2), vertcat(diff(autotime)/2, 0)),2); % in ns
    autowidth = autowidth * 1e-9;   % in s
    
    NCounts = head.Records;
 
    [sync, micro, chan, special] = ht3v2read_old(name, [1, photons]);
    tchan = chan(special == 0);
    dind = unique(tchan(1:1e4));
    ndet = length(dind);
    dindrev_dind = zeros(8,1);
    for j=1:length(dind)
        rev_dind(dind(j)+1) = j;
    end

    auto = zeros(length(autotime), ndet, ndet, 100);
    auto_tmp = zeros(length(autotime), ndet, ndet, 10);

    auto2 = zeros(length(autotime2), ndet, ndet, 100);
    auto2_tmp = zeros(length(autotime2), ndet, ndet, 10);

    clf;
    tcspc = zeros(2^15, ndet);
    rate = zeros(ceil(NCounts/photons/10), ndet);
    rate_tmp = zeros(10, ndet);
    time = zeros(ceil(NCounts/photons/10), 1);
    time_tmp = zeros(10, 1);
    jj = 1;
    jj_tmp = 1;
    jj_binned = 0;
    cnt = 0;
    tst = 1;
    sync2 = [];
    micro2 = [];
    chan2 = [];    
    special2 = [];
    overcount2 = 0;
    while tst>0
        [sync, micro, chan, special, tst, overcount] = ht3v2read_old(name, [cnt+1, photons]);
        if tst>0
            cnt = cnt + tst;
           
            sync2 =  vertcat(sync2, sync+overcount2);
            overcount2 = overcount2 + overcount;
            micro2 = vertcat(micro2, micro);
            chan2 = vertcat(chan2, chan);
            special2 = vertcat(special2, special);

            trigposs = find(chan2==4 & special2==1);
            trigposs = vertcat(0, trigposs);                      
            
            for trigcount=2:length(trigposs)
                from = trigposs(trigcount-1)+1;
                to = trigposs(trigcount)-1;
                sync = sync2(from:to);
                micro = micro2(from:to);
                chan = chan2(from:to);
                special = special2(from:to);
                
                sync(special==1) = [];
                micro(special==1) = [];
                chan(special==1) = [];

                if (length(sync) > 1)

                    for j=1:length(chan)
                        dd = rev_dind(chan(j)+1);
                        tcspc(micro(j)+1, dd) = tcspc(micro(j)+1, dd) + 1;
                    end

                    tim = sync*syncperiod + micro*microresolution;  % photon arrival time in ns
                    time_tmp(jj_tmp) = (tim(end) - tim(1))*1e-9;            % length in s

                    for k=1:ceil(maxtime(1)/75)
                        tmp = tim(k+1:end)-tim(1:end-k);
                        for det1 = 1:ndet
                            for det2 = 1:ndet
                                ind = chan(1:end-k)==dind(det1) & chan(k+1:end)==dind(det2);
                                if sum(ind)>1
                                    auto_tmp(:,det1,det2,jj_tmp) = auto_tmp(:,det1,det2,jj_tmp) + mHist(tmp(ind),autotime)./autowidth;
                                end
                            end
                        end
                    end

                    for det1 = 1:ndet
                        rate_tmp(jj_tmp,det1) = sum(chan==dind(det1))/time_tmp(jj_tmp); % rate per second
                    end

                    flvp = zeros(length(chan),ndet);
                    for p=1:ndet
                        flvp(:,p) = chan==dind(p);
                    end
                    auto2_tmp(:,:,:,jj_tmp) = tttr2xfcs(sync, flvp, Ncasc2, Nsub) / Timeunit2;

                    for det1 = 1:ndet
                        for det2 = 1:ndet
                            auto_tmp(:,det1,det2,jj_tmp) = ((auto_tmp(:,det1,det2,jj_tmp) / time_tmp(jj_tmp) / (rate_tmp(jj_tmp,det1) * rate_tmp(jj_tmp,det2)))) - 1;
                            auto2_tmp(:,det1,det2,jj_tmp) = ((auto2_tmp(:,det1,det2,jj_tmp) / time_tmp(jj_tmp) / (rate_tmp(jj_tmp,det1) * rate_tmp(jj_tmp,det2)))) - 1;
                        end
                    end

                    if (jj_tmp >= 10)
                        jj_tmp = 0;
                        jj_binned = jj_binned + 1;
                        auto(:,:,:,jj_binned) = mean(auto_tmp, 4); 
                        auto2(:,:,:,jj_binned) = mean(auto2_tmp, 4);
                        time(jj_binned) = sum(time_tmp);
                        rate(jj_binned,:) = mean(rate_tmp,1);
                        
                        auto_tmp = zeros(length(autotime), ndet, ndet, 10);
                        auto2_tmp = zeros(length(autotime2), ndet, ndet, 10);
                        rate_tmp = zeros(10, ndet);
                        time_tmp = zeros(10, 1);
                    end

                    subplot('position', [0.925 0.2 0.025 0.6]);
                    bar(0,cnt/NCounts);
                    axis([-0.4 0.4 0 1]);
                    set(gca,'xtick',[],'ytick',[]);
                    subplot('position', [0.1 0.1 0.7 0.8]);

                    if (jj_binned > 1)
                        blaa.auto = auto;
                        blaa.autotime = autotime*1e-9;
                        blaa2.auto = auto2;
                        blaa2.autotime = autotime2;
                        [y, t] = FCSCrossReadCW(blaa);
                        [y2, t2] = FCSCrossReadCW(blaa2);
                        semilogx(t, mean(y(:,:,1:jj_binned),3), t2(3:end), mean(y2(3:end,:,1:jj_binned),3));
                        ylim([0 max(max(mean(y(:,:,1:jj_binned), 3)))*1.05]);
                    end
                    drawnow;
                        
                    jj = jj+1;
                    jj_tmp = jj_tmp + 1;
                    sizes = size(auto);
                    sizes2 = size(auto2);
                    if (jj_binned > sizes(4))
                        auto = cat(4, auto, zeros([sizes(1:3) 100])); 
                        auto2 = cat(4, auto2, zeros([sizes2(1:3) 100])); 
                    end
                end
            end
            fro = trigposs(end)+1;
            sync2 = sync2(fro:end);
            micro2 = micro2(fro:end);
            chan2 = chan2(fro:end);
            special2 = special2(fro:end);      
        end
    end
    clf

    auto = auto(:,:,:,1:jj_binned);
    auto2 = auto2(:,:,:,1:jj_binned);
    blaa.auto = auto;
    blaa.autotime = autotime*1e-9;
    blaa2.auto = auto2;
    blaa2.autotime = autotime2;
    [y, t] = FCSCrossReadCW(blaa);
    [y2, t2] = FCSCrossReadCW(blaa2);
    semilogx(t, mean(y,3), t2(3:end), mean(y2(3:end,:,:),3));
    %semilogx(autotime*1e-9, reshape(sum(auto,4),[size(auto,1) size(auto,2)*size(auto,3)])/sum(time));

    tau = (0:2^15-1)*microresolution;
    tau = tau';
    tau(any(tcspc==0,2),:)=[];
    tcspc(any(tcspc==0,2),:)=[];

    res.tau = tau;
    res.tcspc = tcspc;

    res.autotime1 = autotime;
    res.autotime2 = autotime2;
    res.auto1 = auto;
    res.auto2 = auto2;

    [res.autotime ind] = sort(vertcat(autotime*1e-9, autotime2(3:end)));   % in s
    tmp = cat(1, auto(:,:,:,:), auto2(3:end,:,:,:));
    res.auto = tmp(ind,:,:,:);

    res.autowidth1 = autowidth;
    res.autowidth = sum(horzcat(vertcat(0, diff(res.autotime)/2), vertcat(diff(res.autotime)/2, 0)),2); % in s

    res.rate = rate(1:jj_binned,:);
    res.time = time(1:jj_binned);    
else
    
end



