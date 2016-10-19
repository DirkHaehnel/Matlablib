function [auto, auto2, tcspc, rate, time] = SingleFocusCW2FCS_Cluster_Calc(filename)

    load(filename);
    
    tcspc = zeros(2^15, ndet);
    sync(special==1) = [];
    micro(special==1) = [];
    chan(special==1) = [];
    
    auto = zeros(length(autotime), ndet, ndet);
    auto2 = zeros(length(autotime2), ndet, ndet);
    rate = zeros(ndet,1);

    for j=1:length(chan)
        dd = rev_dind(chan(j)+1);
        tcspc(micro(j)+1, dd) = tcspc(micro(j)+1, dd) + 1;
    end

    tim = sync*syncperiod + micro*microresolution;  % photon arrival time in ns
    time = (tim(end) - tim(1))*1e-9;                % length in s

    for k=1:ceil(maxtime(1)/75)
        tmp = tim(k+1:end)-tim(1:end-k);
        for det1 = 1:ndet
            for det2 = 1:ndet
                ind = chan(1:end-k)==dind(det1) & chan(k+1:end)==dind(det2);
                if sum(ind)>1
                    auto(:,det1,det2) = auto(:,det1,det2) + mHist(tmp(ind),autotime)./autowidth;
                end
            end
        end
    end

    for det1 = 1:ndet
        rate(det1) = sum(chan==dind(det1))/time; % rate per second
    end

    flvp = zeros(length(chan),ndet);
    for p=1:ndet
        flvp(:,p) = chan==dind(p);
    end
    auto2(:,:,:) = tttr2xfcs(sync, flvp, Ncasc2, Nsub) / Timeunit2;

    for det1 = 1:ndet
        for det2 = 1:ndet
            auto(:,det1,det2) = ((auto(:,det1,det2) / time / (rate(det1) * rate(det2)))) - 1;
            auto2(:,det1,det2) = ((auto2(:,det1,det2) / time / (rate(det1) * rate(det2)))) - 1;
        end
    end

end
