function [res, head] = ingopt3Pro(name, taskflag, param, weight)

warning off MATLAB:DeprecatedLogicalAPI

timeframe = 1e10;  % time-frame to read photons from file in ns
smooth    = 1;     % smoothed autocorrelation
z32       = 2^32-1;
close all

k = strfind(name, '.');
if isempty(k) 
    name = [name '.ingopt3'];
end;

head       = ingopt3Read_head(name);

if ~isempty(head)
    CntRate0   = head.CntRate0;
    t_Counts   = head.NCounts;
    Resolution = mean(head.Resolution);
    sync       = 1e9/CntRate0;
    NChannels  = 4096;
    end_time   = 0;
 
    if nargin<4
        weight = [];
    end;
    if nargin<3
        param  = [];
    end;
  
    ratefilter = 0;
    p = strfind(taskflag, '_');
    if numel(p)>0
       taskflag = taskflag(1:p(1)-1);
       ratefilter = 1;
    end


    if strcmp(taskflag,'ratefilter')
        
        if ~isempty(param)
            timeframe = param(1);
        end;
        
        rate    = [];
        time    = [];
        tst     = 1;
        cnt     = 0;
        cnt1    = 0;

        max_time  = 0;
        overcount = 0;

        h = waitbar(0,'Please wait ...');
        while (tst>0)

            max_time = max_time + timeframe;

            [tmpy, tmpx, flag, tst, overcount] = ingopt3Read_t(name, cnt1+1, overcount, max_time, sync);

            if (tst>0)

                cnt1 = cnt1 + tst;

                phot = numel(flag(flag==1|flag==2));
                dt   = (tmpy(end)-tmpy(1)+(tmpx(end)-tmpx(1))*Resolution)*1e-9;
                cnt  = cnt+phot;
                time = [time dt];
                rate = [rate phot/dt];
            end
            waitbar(cnt1/t_Counts);
            drawnow;
        end
        close(h);

        corrind = ones(size(rate));
        if ratefilter
            mrate   = mean(rate);
            corrind = abs(rate-mrate)./mrate<0.05;
        end;
        
        res.rate    = rate;
        res.time    = time;
        res.cnt     = cnt;
        res.corrind = corrind;
    end;
    

    if strcmp(taskflag,'tcspc')
        
        timeunit = sync;
        
        cnt     = 0;
        corrind = [];
        time    = [];
        rate    = [];
        
        if ratefilter
            res     = ingopt3Pro(name, 'ratefilter');
            corrind = res.corrind;
            cnt     = res.cnt;
            time    = res.time;
            rate    = res.rate;
        end
        
        tau     = (0:(NChannels-1));
        tcspc   = zeros(NChannels,2);
            
        tst       = 1;
        cnt1      = 0;
        max_time  = 0;
        overcount = 0;
        n         = 1;

        while (tst>0)

            max_time = max_time + timeframe;

            [tmpy, tmpx, flag, tst, overcount] = ingopt3Read_t(name, cnt1+1, overcount, max_time, sync);

            if (tst>0) 

                cnt1 = cnt1 + tst;
                
                if numel(corrind)<n
                    corrind = [corrind 1];
                    phot = numel(flag(flag==1|flag==2));
                    dt   = (tmpy(end)-tmpy(1)+(tmpx(end)-tmpx(1))*Resolution)*1e-9;
                    cnt  = cnt+phot;
                    time = [time dt];
                    rate = [rate phot/dt];
                end;
                
                if corrind(n)                                        
                    tcspc(:,1) = tcspc(:,1) + mHist(tmpx(flag==1), tau);
                    tcspc(:,2) = tcspc(:,2) + mHist(tmpx(flag==2), tau);

                    subplot('position', [0.875 0.2 0.05 0.6]);
                    bar(0, cnt1/t_Counts);
                    axis([-0.4 0.4 0 1]);
                    set(gca,'xtick',[],'ytick',[]);
                    subplot('position', [0.1 0.1 0.7 0.8]);
                    semilogy(tau*Resolution,tcspc);
                    drawnow;
                end;
                n = n + 1;
            end

        end

        semilogy(tau*Resolution, tcspc);

        res.tau   = tau(sum(tcspc,2)>0)*Resolution;
        res.tcspc = tcspc(sum(tcspc,2)>0,:);
    end


    if strcmp(taskflag,'bin')

        timeframe = 1e9;  % time-frame to read photons from file in ns
        timeunit  = sync;

        binwidth = 1e6/timeunit;      % in timeunits
        if ~isemty(param)      
            binwidth = param(1)/timeunit;  % in timeunits
        end

        rate1   = [];
        rate2   = [];
        time    = [];
        bin1    = [];
        bin2    = [];
        cnt     = 0;
        cnt1    = 0;
        tst     = 1;

        max_time  = 0;
        overcount = 0;

        while (tst>0)

            max_time = max_time + timeframe;

            [tmpy, tmpx, flag, tst, overcount] = ingopt3Read_t(name, cnt1+1, overcount, max_time, sync);

            if (tst>0)

                cnt1 = cnt1 + tst;

                tmpy = round((tmpy+tmpx*Resolution)/timeunit);  % in timeunits
                tmpy = tmpy-tmpy(1);

                dt    = tmpy(end)*timeunit*1e-9;
                time  = [time dt];

                rate1 = [rate1 length(tmpy(flag==1))/dt];
                rate2 = [rate2 length(tmpy(flag==2))/dt];

                tmp1 = tttr2bin(tmpy(flag==1), binwidth);
                tmp2 = tttr2bin(tmpy(flag==2), binwidth);
                bin1 = [bin1 tmp1'];
                bin2 = [bin2 tmp2'];

                ind = (1:min([length(tmp1) length(tmp2)]));

                subplot('position', [0.875 0.2 0.05 0.6]);
                bar(0,cnt1/t_Counts);
                axis([-0.4 0.4 0 1]);
                set(gca,'xtick',[],'ytick',[]);
                subplot('position', [0.1 0.1 0.7 0.8]);
                bar(1e-9*ind*binwidth*sync,tmp2(ind)+tmp1(ind));
                axis tight
                drawnow;
            end

        end
        clf

        res.bin1  = bin1;
        res.bin2  = bin2;
        res.rate1 = rate1;
        res.rate2 = rate2;
        res.time  = time;
    end


    if strcmp(taskflag,'fcs')   % autocorrelation of each detector channel

        
        Nsub     = 10;
        timeunit = sync;
        mct      = 3;          % 3 sec max correlation time
        if ~isempty(param)
            mct   = param(1);
            if length(param)>1
                Nsub  = param(2);
            end;
            if length(param)>2
                timeunit = param(3)*Resolution;
            end
        end

        Ncasc    = ceil(log2(mct/(timeunit*1e-9*Nsub)+1));
        autotime = 1e-9*timeunit*cumsum(reshape(repmat(2.^[0 0:Ncasc-1],Nsub,1),(Ncasc+1)*Nsub,1));
        auto     = zeros(length(autotime),3);

        tau   = (0:(NChannels-1));
        tcspc = zeros(NChannels,2);

        rate    = [];
        time    = [];
        tst     = 1;
        cnt     = 0;
        cnt1    = 0;

        max_time  = 0;
        overcount = 0;

        while (tst>0)

            max_time = max_time + timeframe;
            if end_time > 0
                if max_time > end_time
                    max_time = end_time;
                end;
            end;

            [tmpy, tmpx, flag, tst, overcount] = ingopt3Read_t(name, cnt1+1, overcount, max_time, sync);

            if (tst>0)

                cnt1 = cnt1 + tst;

                tmpy = round((tmpy+tmpx*Resolution)/timeunit);  % in timeunits
                tmpy = tmpy-tmpy(1);

                cnt = cnt+length(tmpy);
                dt =  tmpy(end)*timeunit*1e-9;  % in s
                rate = [rate length(tmpy)/dt];
                time = [time dt];

                tcspc(:,1) = tcspc(:,1) + mHist(tmpx(flag==1), tau);
                tcspc(:,2) = tcspc(:,2) + mHist(tmpx(flag==2), tau);

                if tmpy(end)>z32
                    auto(:,1) = auto(:,1) + ingopt32fcs(tmpy, double(flag==1), double(flag==1), Ncasc, Nsub, smooth);
                    auto(:,2) = auto(:,2) + ingopt32fcs(tmpy, double(flag==2), double(flag==2), Ncasc, Nsub, smooth);
                    auto(:,3) = auto(:,3) + ingopt32fcs(tmpy, double(flag<=2), double(flag<=2), Ncasc, Nsub, smooth);
                else
                    auto(:,1) = auto(:,1) + tttr2fcs(tmpy, double(flag==1), double(flag==1), Ncasc, Nsub, smooth);
                    auto(:,2) = auto(:,2) + tttr2fcs(tmpy, double(flag==2), double(flag==2), Ncasc, Nsub, smooth);
                    auto(:,3) = auto(:,3) + tttr2fcs(tmpy, double(flag<=2), double(flag<=2), Ncasc, Nsub, smooth);
                end

                subplot('position', [0.875 0.2 0.05 0.6]);
                bar(0,cnt1/t_Counts);
                axis([-0.4 0.4 0 1]);
                set(gca,'xtick',[],'ytick',[]);
                subplot('position', [0.1 0.1 0.7 0.8]);
                semilogx(autotime,auto/sum(time)/timeunit*1e9);
                drawnow;
            end

        end
        clf

        semilogx(autotime,auto/sum(time)/timeunit*1e9);

        res.autotime = autotime;
        res.auto     = auto/sum(time)/timeunit*1e9;
        res.tau      = tau(sum(tcspc,2)>0)*Resolution;
        res.tcspc    = tcspc(sum(tcspc,2)>0,:);
        res.rate     = rate;
        res.time     = time;

    end

    if strcmp(taskflag,'fcs3')   % autocorrelation of each detector channel

        
        Nsub     = 10;
        timeunit = sync;
        mct      = 3;          % 3 sec max correlation time
        if ~isempty(param)
            mct   = param(1);
            if length(param)>1
                Nsub  = param(2);
            end;
            if length(param)>2
                timeunit = param(3)*Resolution;
            end
        end

        Ncasc    = ceil(log2(mct/(timeunit*1e-9*Nsub)+1));
        autotime = 1e-9*timeunit*cumsum(reshape(repmat(2.^[0 0:Ncasc-1],Nsub,1),(Ncasc+1)*Nsub,1));
        auto     = zeros(length(autotime),1);

        tau   = (0:(NChannels-1));
        tcspc = zeros(NChannels,1);

        rate    = [];
        time    = [];
        tst     = 1;
        cnt     = 0;
        cnt1    = 0;

        max_time  = 0;
        overcount = 0;

        while (tst>0)

            max_time = max_time + timeframe;
            if end_time > 0
                if max_time > end_time
                    max_time = end_time;
                end;
            end;

            [tmpy, tmpx, flag, tst, overcount] = ingopt3Read_t(name, cnt1+1, overcount, max_time, sync);

            if (tst>0)

                cnt1 = cnt1 + tst;

                tmpy = round((tmpy+tmpx*Resolution)/timeunit);  % in timeunits
                tmpy = tmpy-tmpy(1);

                cnt = cnt+length(tmpy);
                dt =  tmpy(end)*timeunit*1e-9;  % in s
                rate = [rate length(tmpy)/dt];
                time = [time dt];

                tcspc = tcspc + mHist(tmpx(flag==3), tau);

                if tmpy(end)>z32
                    auto(:,1) = auto(:,1) + ingopt32fcs(tmpy, double(flag==3), double(flag==3), Ncasc, Nsub, smooth);
                else
                    auto(:,1) = auto(:,1) + tttr2fcs(tmpy, double(flag==3), double(flag==3), Ncasc, Nsub, smooth);
                end

                subplot('position', [0.875 0.2 0.05 0.6]);
                bar(0,cnt1/t_Counts);
                axis([-0.4 0.4 0 1]);
                set(gca,'xtick',[],'ytick',[]);
                subplot('position', [0.1 0.1 0.7 0.8]);
                semilogx(autotime,auto/sum(time)/timeunit*1e9);
                drawnow;
            end

        end
        clf

        semilogx(autotime,auto/sum(time)/timeunit*1e9);

        res.autotime = autotime;
        res.auto     = auto/sum(time)/timeunit*1e9;
        res.tau      = tau(sum(tcspc,2)>0)*Resolution;
        res.tcspc    = tcspc(sum(tcspc,2)>0,:);
        res.rate     = rate;
        res.time     = time;

    end

    
    if strcmp(taskflag,'xfcs')
        
        Nsub     = 10;
        timeunit = sync;
        mct      = 3;          % 3 sec max correlation time
        filter   = [1 sync/Resolution];
        
        if ~isempty(param)
            mct   = param(1);
            if length(param)>1
                Nsub  = param(2);
            end;
            if length(param)>2
                timeunit = param(3)*Resolution;
            end
        end

        if ~isempty(weight)
            filter = weight./Resolution;
        end
        
        cnt     = 0;
        corrind = [];
        time    = [];
        rate    = [];
        
        if ratefilter
            res     = ingopt3Pro(name, 'ratefilter');
            corrind = res.corrind;
            cnt     = res.cnt;
            time    = res.time;
            rate    = res.rate;
        end
        
        Ncasc    = ceil(log2(mct/(timeunit*1e-9*Nsub)+1));
        autotime = 1e-9*timeunit*cumsum(reshape(repmat(2.^[0 0:Ncasc-1],Nsub,1),(Ncasc+1)*Nsub,1));
        auto     = zeros(length(autotime),2);

        tau   = (0:(NChannels-1));
        tcspc = zeros(NChannels,2);
        
        tst       = 1;
        cnt1      = 0;
        max_time  = 0;
        overcount = 0;
        n         = 1;

        while (tst>0)

            max_time = max_time + timeframe;

            [tmpy, tmpx, flag, tst, overcount] = ingopt3Read_t(name, cnt1+1, overcount, max_time, sync);

            if (tst>0) 

                cnt1 = cnt1 + tst;
                
                if numel(corrind)<n
                    corrind = [corrind 1];
                    phot = numel(flag(flag==1|flag==2));
                    dt   = (tmpy(end)-tmpy(1)+(tmpx(end)-tmpx(1))*Resolution)*1e-9;
                    cnt  = cnt+phot;
                    time = [time dt];
                    rate = [rate phot/dt];
                end;

                if corrind(n)                                        

                    ind = (tmpx<filter(1)|tmpx>filter(2));

                    tmpy(ind) = [];
                    tmpx(ind) = [];
                    flag(ind) = [];

                    tmpy = round((tmpy+tmpx*Resolution)/timeunit);  % in timeunits
                    tmpy = tmpy-tmpy(1);

                    tcspc(:,1) = tcspc(:,1) + mHist(tmpx(flag==1), tau);
                    tcspc(:,2) = tcspc(:,2) + mHist(tmpx(flag==2), tau);

                    if tmpy(end)>z32
                        auto(:,1) = auto(:,1) + ingopt32fcs(double(tmpy), double(flag==1), double(flag==2), Ncasc, Nsub, smooth);
                        auto(:,2) = auto(:,2) + ingopt32fcs(double(tmpy), double(flag==2), double(flag==1), Ncasc, Nsub, smooth);
                    else
                        auto(:,1) = auto(:,1) + tttr2fcs(double(tmpy), double(flag==1), double(flag==2), Ncasc, Nsub, smooth);
                        auto(:,2) = auto(:,2) + tttr2fcs(double(tmpy), double(flag==2), double(flag==1), Ncasc, Nsub, smooth);
                    end

                    subplot('position', [0.875 0.2 0.05 0.6]);
                    bar(0, cnt1/t_Counts);
                    axis([-0.4 0.4 0 1]);
                    set(gca,'xtick',[],'ytick',[]);
                    subplot('position', [0.1 0.1 0.7 0.8]);
                    semilogx(autotime,auto/sum(time(corrind==1))/timeunit*1e9);
                    drawnow;
                end;
                n = n + 1;
            end

        end
        clf

        semilogx(autotime,auto/sum(time(corrind==1))/timeunit*1e9);

        res.autotime = autotime;
        res.auto     = auto/sum(time(corrind==1))/timeunit*1e9;
        res.tau      = tau(sum(tcspc,2)>0)*Resolution;
        res.tcspc    = tcspc(sum(tcspc,2)>0,:);
        res.rate     = rate;
        res.time     = time;
    end


    if strcmp(taskflag,'xfcsm')         % cross-correlation of individual time-frames
        
        Nsub     = 10;
        timeunit = sync;
        mct      = 3;          % 3 sec max correlation time
        filter   = [1 sync/Resolution];
        
        if ~isempty(param)
            mct   = param(1);
            if length(param)>1
                Nsub  = param(2);
            end;
            if length(param)>2
                timeunit = param(3)*Resolution;
            end
        end

        if ~isempty(weight)
            filter = weight./Resolution;
        end

        tau   = (0:(NChannels-1));
        tcspc = zeros(NChannels,2);

        rate   = [];
        time   = [];
        tst    = 1;
        cnt    = 0;
        cnt1   = 0;
        scnt   = 0;

        max_time  = 0;
        overcount = 0;

        while (tst>0)

            max_time = max_time + timeframe;

            [tmpy, tmpx, flag, tst, overcount] = ingopt3Read_t(name, cnt1+1, overcount, max_time, sync);


            if (tst>0)

                cnt1 = cnt1 + tst;
                scnt = scnt + 1;

                phot = numel(flag(flag==1|flag==2))

                ind       = (tmpx<filter(1)|tmpx>filter(2));
                tmpy(ind) = [];
                tmpx(ind) = [];
                flag(ind) = [];

                tmpy = round((tmpy+tmpx*Resolution)/timeunit);  % in timeunits
                tmpy = tmpy-tmpy(1);

                cnt  = cnt+phot;
                dt   = (tmpy(end))*timeunit*1e-9;
                rate = [rate phot/dt];
                time = [time dt];

                tcspc(:,1,scnt) = mHist(tmpx(flag==1), tau);
                tcspc(:,2,scnt) = mHist(tmpx(flag==2), tau);

                if tmpy(end)>z32
                    auto(:,1,scnt) = ingopt32fcs(tmpy, double(flag==1), double(flag==2), Ncasc, Nsub, smooth);
                    auto(:,2,scnt) = ingopt32fcs(tmpy, double(flag==2), double(flag==1), Ncasc, Nsub, smooth);
                else
                    auto(:,1,scnt) = tttr2fcs(tmpy, double(flag==1), double(flag==2), Ncasc, Nsub, smooth);
                    auto(:,2,scnt) = tttr2fcs(tmpy, double(flag==2), double(flag==1), Ncasc, Nsub, smooth);
                end

                subplot('position', [0.875 0.2 0.05 0.6]);
                bar(0, cnt1/t_Counts);
                axis([-0.4 0.4 0 1]);
                set(gca,'xtick',[],'ytick',[]);
                subplot('position', [0.1 0.1 0.7 0.8]);
                semilogx(autotime,squeeze(sum(auto,3))/sum(time)/timeunit*1e3);
                drawnow;
            end

        end
        clf

        semilogx(autotime,squeeze(sum(auto,3))/sum(time)/timeunit*1e3);
        res.autotime = autotime;
        res.auto     = auto/sum(time)/timeunit*1e3;
        res.tau      = tau(sum(tcspc,2)>0)*Resolution;
        res.tcspc    = tcspc(sum(tcspc,2)>0,:);
        res.rate     = rate;
        res.time     = time;

    end


    if strcmp(taskflag,'xfcslin')    % Cross-Correlation with linear timescale
        
        Nsub     = 500;
        timeunit = sync;
        filter   = [1 sync/Resolution];
        
        if ~isempty(param)
            Nsub   = param(1);
            if length(param)>1
                timeunit = param(2)*Resolution;
            end;
        end

        if ~isempty(weight)
            filter = weight./Resolution;
        end

        cnt     = 0;
        corrind = [];
        time    = [];
        rate    = [];
        
        if ratefilter
            res     = ingopt3Pro(name, 'ratefilter');
            corrind = res.corrind;
            cnt     = res.cnt;
            time    = res.time;
            rate    = res.rate;
        end

        Ncasc    = 1;
        autotime = 1e-9*timeunit*(1:1000)';
        auto     = zeros(length(autotime),2);

        tau   = (0:(NChannels-1));
        tcspc = zeros(NChannels,2);

        tst       = 1;
        cnt1      = 0;
        max_time  = 0;
        overcount = 0;
        n         = 1;

        while (tst>0)

            max_time = max_time + timeframe;

            [tmpy, tmpx, flag, tst, overcount] = ingopt3Read_t(name, cnt1+1, overcount, max_time, sync);

            if (tst>0) 

                cnt1 = cnt1 + tst;

                if numel(corrind)<n
                    corrind = [corrind 1];
                    phot = numel(flag(flag==1|flag==2));
                    dt   = (tmpy(end)-tmpy(1)+(tmpx(end)-tmpx(1))*Resolution)*1e-9;
                    cnt  = cnt+phot;
                    time = [time dt];
                    rate = [rate phot/dt];
                end;

                if corrind(n)                                        
                    ind = (tmpx<filter(1)|tmpx>filter(2));

                    tmpy(ind) = [];
                    tmpx(ind) = [];
                    flag(ind) = [];

                    tmpy = round((tmpy+tmpx*Resolution)/timeunit);  % in timeunits
                    tmpy = tmpy-tmpy(1);

                    tcspc(:,1) = tcspc(:,1) + mHist(tmpx(flag==1), tau);
                    tcspc(:,2) = tcspc(:,2) + mHist(tmpx(flag==2), tau);

                    if tmpy(end)>z32
                        auto(:,1) = auto(:,1) + ingopt32fcs(double(tmpy), double(flag==1), double(flag==2), Ncasc, Nsub, smooth);
                        auto(:,2) = auto(:,2) + ingopt32fcs(double(tmpy), double(flag==2), double(flag==1), Ncasc, Nsub, smooth);
                    else
                        auto(:,1) = auto(:,1) + tttr2fcs(double(tmpy), double(flag==1), double(flag==2), Ncasc, Nsub, smooth);
                        auto(:,2) = auto(:,2) + tttr2fcs(double(tmpy), double(flag==2), double(flag==1), Ncasc, Nsub, smooth);
                    end

                    subplot('position', [0.875 0.2 0.05 0.6]);
                    bar(0, cnt1/t_Counts);
                    axis([-0.4 0.4 0 1]);
                    set(gca,'xtick',[],'ytick',[]);
                    subplot('position', [0.1 0.1 0.7 0.8]);
                    semilogx(autotime,auto/sum(time(corrind==1))/timeunit*1e9);
                    drawnow;
                end;
                n = n + 1;
                
            end

        end
        clf

        auto = auto/sum(time(corrind==1))/timeunit*1e9;

        plot(autotime, mean(auto,2));

        res.autotime = autotime;
        res.auto     = auto;
        res.tau      = tau(sum(tcspc,2)>0)*Resolution;
        res.tcspc    = tcspc(sum(tcspc,2)>0,:);
        res.rate     = rate;
        res.time     = time;

    end


    if strcmp(taskflag,'2ffcs')      % 2-focus FCS

        Nsub     = 10;
        timeunit = sync;
        mct      = 3;          % 3 sec max correlation time
        filter   = [1 sync/Resolution];
        
        if ~isempty(param)
            mct   = param(1);
            if length(param)>1
                Nsub  = param(2);
            end;
            if length(param)>2
                timeunit = param(3)*Resolution;
            end
        end

        if ~isempty(weight)
            filter = weight./Resolution;
        end

        cnt     = 0;
        corrind = [];
        time    = [];
        rate    = [];
        
        if ratefilter
            res     = ingopt3Pro(name, 'ratefilter');
            corrind = res.corrind;
            cnt     = res.cnt;
            time    = res.time;
            rate    = res.rate;
        end

        
        Ncasc     = ceil(log2(mct/(timeunit*1e-9*Nsub)+1));
        autotime  = 1e-9*timeunit*cumsum(reshape(repmat(2.^[0 0:Ncasc-1],Nsub,1),(Ncasc+1)*Nsub,1));
        auto      = zeros(length(autotime), 4);  
        tau       = (0:(NChannels-1));
        tcspc     = zeros(NChannels, 4);
        tcspcdata = zeros(NChannels, 2);

        tst       = 1;
        cnt1      = 0;
        overcount = 0;

        while (tst>0)

            [tmpy, tmpx, flag, tst, overcount] = ingopt3Read(name, cnt1+1, 5e5, overcount, sync);

            if (tst>0)
                cnt1 = cnt1 + tst;

                tcspcdata(:,1) = tcspcdata(:,1) + mHist(tmpx(flag==1), tau);
                tcspcdata(:,2) = tcspcdata(:,2) + mHist(tmpx(flag==2), tau);

                subplot('position', [0.925 0.2 0.025 0.6]);
                bar(0, cnt1/t_Counts);
                axis([-0.4 0.4 0 1]);
                set(gca,'xtick',[],'ytick',[]);
                subplot('position', [0.1 0.1 0.7 0.8]);
                semilogy(tau*Resolution,tcspcdata);
                drawnow;
            end

        end

        ind = sum(tcspcdata,2)==0;
        tcspcdata(ind,:) = [];

        win = floor(sync/Resolution/2);      % Max #channels / pulse;

        for d = 1:2                          % detectors

            [pos pos1] = max(tcspcdata(    1:win, d));
            [pos pos2] = max(tcspcdata(win+1:end, d));
            pos2 = pos2 + win;
            x1   = pos1 - 1;
            x2   = size(tcspcdata, 1) - pos2 - 1;

            if ((pos1 + x2)>(pos2 - x1))
                if (x1>=x2)
                    x1 = pos2-pos1-x1-1;
                else
                    x2 = pos2-pos1-x1-1;
                end;
            end;

            tind(d, 1, 1) = pos1-x1;
            tind(d, 1, 2) = pos1+x2;
            tind(d, 2, 1) = pos2-x1;
            tind(d, 2, 2) = pos2+x2;

            nanz = 1+x2+x1;
            tcspc(1:nanz,2*d-1) = tcspcdata((pos1-x1:pos1+x2), d);
            tcspc(1:nanz,2*d)   = tcspcdata((pos2-x1:pos2+x2), d);

        end

        ind = sum(tcspc,2)==0;
        tcspc(ind,:) = [];
        tau = 0:size(tcspc,1)-1;
        semilogy(tau*Resolution, tcspc); drawnow

        res.tau       = Resolution*tau;
        res.tcspc     = tcspc;
        res.tcspcdata = tcspcdata;

        tst       = 1;
        cnt1      = 0;
        max_time  = 0;
        overcount = 0;
        n         = 1;

        while (tst>0)
            max_time = max_time + timeframe;

            [tmpy, tmpx, flag, tst, overcount] = ingopt3Read_t(name, cnt1+1, overcount, max_time, sync);

            if (tst>0)

                cnt1 = cnt1 + tst;

                if numel(corrind)<n
                    corrind = [corrind 1];
                    phot = numel(flag(flag==1|flag==2));
                    dt   = (tmpy(end)-tmpy(1)+(tmpx(end)-tmpx(1))*Resolution)*1e-9;
                    cnt  = cnt+phot;
                    time = [time dt];
                    rate = [rate phot/dt];
                end;

                if corrind(n)                                        
                    ind = (tmpx<filter(1)|tmpx>filter(2));

                    tmpy(ind) = [];
                    tmpx(ind) = [];
                    flag(ind) = [];

                    tmpy = round((tmpy+tmpx*Resolution)/timeunit);  % in timeunits
                    tmpy = tmpy-tmpy(1);

                    fla = double(flag==1 & tmpx>=tind(1, 1, 1) & tind(1, 1, 2)>=tmpx);
                    flb = double(flag==2 & tmpx>=tind(2, 1, 1) & tind(2, 1, 2)>=tmpx);
                    flc = double(flag==2 & tmpx>=tind(2, 1, 1) & tind(2, 1, 2)>=tmpx);
                    fld = double(flag==1 & tmpx>=tind(1, 1, 1) & tind(1, 1, 2)>=tmpx);

                    auto(:, 1, n) = (tttr2fcs(tmpy, fla, flb, Ncasc, Nsub, smooth)+tttr2fcs(tmpy, flc, fld, Ncasc, Nsub, smooth))/dt/timeunit;

                    fla = double(flag==1 & tmpx>=tind(1, 1, 1) & tind(1, 1, 2)>=tmpx);
                    flb = double(flag==2 & tmpx>=tind(2, 2, 1) & tind(2, 2, 2)>=tmpx);
                    flc = double(flag==2 & tmpx>=tind(2, 2, 1) & tind(2, 2, 2)>=tmpx);
                    fld = double(flag==1 & tmpx>=tind(1, 1, 1) & tind(1, 1, 2)>=tmpx);

                    auto(:, 2, n) = (tttr2fcs(tmpy, fla, flb, Ncasc, Nsub, smooth)+tttr2fcs(tmpy, flc, fld, Ncasc, Nsub, smooth))/dt/timeunit;

                    fla = double(flag==1 & tmpx>=tind(1, 2, 1) & tind(1, 2, 2)>=tmpx);
                    flb = double(flag==2 & tmpx>=tind(2, 1, 1) & tind(2, 1, 2)>=tmpx);
                    flc = double(flag==2 & tmpx>=tind(2, 1, 1) & tind(2, 1, 2)>=tmpx);
                    fld = double(flag==1 & tmpx>=tind(1, 2, 1) & tind(1, 2, 2)>=tmpx);

                    auto(:, 3, n) = (tttr2fcs(tmpy, fla, flb, Ncasc, Nsub, smooth)+tttr2fcs(tmpy, flc, fld, Ncasc, Nsub, smooth))/dt/timeunit;

                    fla = double(flag==1 & tmpx>=tind(1, 2, 1) & tind(1, 2, 2)>=tmpx);
                    flb = double(flag==2 & tmpx>=tind(2, 2, 1) & tind(2, 2, 2)>=tmpx);
                    flc = double(flag==2 & tmpx>=tind(2, 2, 1) & tind(2, 2, 2)>=tmpx);
                    fld = double(flag==1 & tmpx>=tind(1, 2, 1) & tind(1, 2, 2)>=tmpx);

                    auto(:, 4, n) = (tttr2fcs(tmpy, fla, flb, Ncasc, Nsub, smooth)+tttr2fcs(tmpy, flc, fld, Ncasc, Nsub, smooth))/dt/timeunit;

                    subplot('position', [0.925 0.2 0.025 0.6]);
                    bar(0,cnt/t_Counts);
                    axis([-0.4 0.4 0 1]);
                    set(gca,'xtick',[],'ytick',[]);
                    subplot('position', [0.1 0.1 0.7 0.8]);
                    semilogx(autotime, [sum(auto(:,1,:),3) sum(auto(:,2,:),3) sum(auto(:,3,:),3) sum(auto(:,4,:),3)]);
                    drawnow;
                end;
                n = n + 1;
            end
        end
        clf
        semilogx(autotime, [sum(auto(:,1,:),3) sum(auto(:,2,:),3) sum(auto(:,3,:),3) sum(auto(:,4,:),3)]);

        drawnow;

        res.autotime = autotime;
        res.auto     = auto;
        res.rate     = rate;
        res.time     = time;

        head.NCounts = cnt;

    end;



    if strcmp(taskflag,'flcs')       % fluorescence lifetime correlation

        if isempty(param)
            disp('No master patterns given!');
            return
        end

        tau   = param(:,1);
        param = param(:,2:end);

        cnt     = 0;
        corrind = [];
        time    = [];
        rate    = [];

        if ratefilter
            if isempty(weight)
                res = ingopt3Pro(name, 'tcspc_');
                weight = res.tcspc(res.tau>=tau(1) & res.tau<=tau(end),:);
            else
                res     = ingopt3Pro(name, 'ratefilter');
            end;
            corrind = res.corrind;
            cnt     = res.cnt;
            time    = res.time;
            rate    = res.rate;
        end

        Nsub  = 10;
        Ncasc = ceil(log2(3/(sync*1e-9*Nsub)+1));    % 3 sec max correlation time
        autotime = 1e-9*sync*cumsum(reshape(repmat(2.^[0 0:Ncasc-1],Nsub,1),(Ncasc+1)*Nsub,1));
        timeunit = sync;
                    
        tst     = 1;
        n       = 1;
        cnt1    = 0;
               
        max_time  = 0;
        overcount = 0;

        % filter generation
        % para = [dye1 channel#1; dye1 channel#2; dye2 channel#1; dye2 channel#2];

        for j=1:size(param,2)
            param(:,j) = param(:,j)/sum(param(:,j));  % normalization
        end

        master = (param(:,1:2:end)'*diag(1./weight(:,1))*param(:,1:2:end))\(diag(1./weight(:,1))*param(:,1:2:end))';
        master = [master; (param(:,2:2:end)'*diag(1./weight(:,2))*param(:,2:2:end))\(diag(1./weight(:,2))*param(:,2:2:end))']';

        auto = zeros(length(autotime),size(master,2)/2,size(master,2)/2);

        while (tst>0)

            max_time = max_time + timeframe;

            [tmpy, tmpx, flag, tst, overcount] = ingopt3Read_t(name, cnt+1, overcount, max_time, sync);

            if (tst>0)

                cnt1 = cnt1 + tst;

                if numel(corrind)<n
                    corrind = [corrind 1];
                    phot = numel(flag(flag==1|flag==2));
                    dt   = (tmpy(end)-tmpy(1)+(tmpx(end)-tmpx(1))*Resolution)*1e-9;
                    cnt  = cnt+phot;
                    time = [time dt];
                    rate = [rate phot/dt];
                end;

                if corrind(n)                                        

                    tmpx = tmpx*Resolution;              % in ns
                    tmpy = round((tmpy+tmpx)/timeunit);  % in timeunits

                    tst = tmpx>=tau(1) & tmpx<=tau(end);
                    tmpy = tmpy(tst);
                    tmpx = tmpx(tst);
                    tmpy = tmpy-tmpy(1);
                    tmpx = round((tmpx-tau(1))/Resolution)+1; % TCSPC-Channel


                    for j=1:size(master,2)/2
                        for k=1:size(master,2)/2
                            if tmpy(end)>z32
                                auto(:,j,k) = auto(:,j,k) + ingopt32fcs(tmpy, master(tmpx,j).*double(flag==1), master(tmpx,end/2+k).*double(flag==2), Ncasc, Nsub, smooth);
                                auto(:,j,k) = auto(:,j,k) + ingopt32fcs(tmpy, master(tmpx,end/2+j).*double(flag==2), master(tmpx,k).*double(flag==1), Ncasc, Nsub, smooth);
                            else
                                auto(:,j,k) = auto(:,j,k) + tttr2fcs(tmpy, master(tmpx,j).*double(flag==1), master(tmpx,end/2+k).*double(flag==2), Ncasc, Nsub, smooth);
                                auto(:,j,k) = auto(:,j,k) + tttr2fcs(tmpy, master(tmpx,end/2+j).*double(flag==2), master(tmpx,k).*double(flag==1), Ncasc, Nsub, smooth);
                            end

                        end
                    end

                    subplot('position', [0.925 0.2 0.025 0.6]);
                    bar(0, cnt1/t_Counts);
                    axis([-0.4 0.4 0 1]);
                    set(gca,'xtick',[],'ytick',[]);
                    subplot('position', [0.1 0.1 0.7 0.8]);
                    semilogx(autotime, auto(:,:,1)/sum(time(corrind==1))/timeunit*1e3);
                    drawnow;
                end;
                n = n + 1;    
            end
        end
        clf

        semilogx(autotime,auto(:,:,1)/sum(time(corrind==1))/timeunit*1e3);

        res.autotime = autotime;
        res.auto     = auto/sum(time(corrind==1))/timeunit*1e3;
        res.tau      = tau;
        res.weight   = weight;
        res.rate     = rate;
        res.master   = master;
        res.time     = time;
    end


    if strcmp(taskflag,'antibunch')
        
        timeunit = Resolution*6.25;
        autotime = (0:2250+1)';
        filter   = [1 sync/Resolution];
        
        if ~isempty(param)
            timeunit = Resolution*param(1);
            autotime = (0:param(2)+1)';
        end
        if ~isempty(weight)
            filter = weight./Resolution;
        end

        cnt     = 0;
        corrind = [];
        time    = [];
        rate    = [];
        
        if ratefilter
            res     = ingopt3Pro(name, 'ratefilter');
            corrind = res.corrind;
            cnt     = res.cnt;
            time    = res.time;
            rate    = res.rate;
        end

        t_shift = timeframe; 

        auto    = zeros(length(autotime),4);
        tau     = (0:(NChannels-1));
        tcspc   = zeros(NChannels,2);

        tst1    = 1;
        cnt1    = 0;
        cnt2    = 0;

        overcount1 = 0;
        overcount2 = 0;

        max_time1 = 0;
        max_time2 = t_shift;
        
        n         = 2;
        if isempty(corrind)
              corrind = [1];
        end;

        [tmpy2, tmpx2, flag2, tst2, overcount2] = ingopt3Read_t(name, cnt2+1, overcount2, max_time2, sync);

        cnt2 = cnt2 - 1;

        while (tst1>0)

            max_time1 = max_time1 + timeframe;
            max_time2 = max_time2 + timeframe;

            [tmpy1, tmpx1, flag1, tst1, overcount1] = ingopt3Read_t(name, cnt1+1, overcount1, max_time1, sync);

            if (tst1>0)

                cnt1 = cnt1 + tst1;
                
                if numel(corrind)<n
                    corrind = [corrind 1];
                    phot = numel(flag1(flag1==1|flag1==2));
                    dt   = (tmpy1(end)-tmpy1(1)+(tmpx1(end)-tmpx1(1))*Resolution)*1e-9;
                    cnt  = cnt+phot;
                    time = [time dt];
                    rate = [rate phot/dt];
                end;

                if corrind(n-1) & corrind(n)
                    
                    [tmpy2, tmpx2, flag2, tst2, overcount2] = ingopt3Read_t(name, cnt2+1, overcount2, max_time2, sync);
        
                    if (tst2>0)
                    
                        cnt2 = cnt2 + tst2;

                        ind1 = (tmpx1<filter(1)|tmpx1>filter(2));
                        ind2 = (tmpx2<filter(1)|tmpx2>filter(2));

                        tmpy1(ind1) = [];
                        tmpx1(ind1) = [];
                        flag1(ind1) = [];
                        tmpy2(ind2) = [];
                        tmpx2(ind2) = [];
                        flag2(ind2) = [];

                        tmpy  = round((tmpy1+tmpx1*Resolution)/timeunit);                                   % in timeunits
                        tmpz1 = round([tmpy1(flag1==1)+tmpx1(flag1==1)*Resolution;
                            tmpy2(flag2==2)+tmpx2(flag2==2)*Resolution-t_shift ]/timeunit);      % in timeunits
                        tmpz2 = round([tmpy2(flag2==1)+tmpx2(flag2==1)*Resolution-t_shift;
                            tmpy1(flag1==2)+tmpx1(flag1==2)*Resolution         ]/timeunit);      % in timeunits

                        ind1  = [zeros(length(tmpy1(flag1==1)),1); ones(length(tmpy2(flag2==2)),1)];
                        ind2  = [zeros(length(tmpy2(flag2==1)),1); ones(length(tmpy1(flag1==2)),1)];

                        [tmpz1, ord] = sort(tmpz1);
                        ind1 = ind1(ord);
                        [tmpz2, ord] = sort(tmpz2);
                        ind2 = ind2(ord);

                        tmpy  = diff(tmpy);
                        tmpz1 = diff(tmpz1);
                        tmpz2 = diff(tmpz2);

                        tcspc(:,1) = tcspc(:,1) + mHist(tmpx1(flag1==1), tau);
                        tcspc(:,2) = tcspc(:,2) + mHist(tmpx1(flag1==2), tau);

                        auto(:,1) = auto(:,1) + mHist( tmpy, autotime, double(flag1(1:end-1)==1).*double(flag1(2:end)==2));
                        auto(:,2) = auto(:,2) + mHist( tmpy, autotime, double(flag1(1:end-1)==2).*double(flag1(2:end)==1));

                        auto(:,3) = auto(:,3) + mHist(tmpz1, autotime, double(ind1(1:end-1)==0).*double(ind1(2:end)==1));
                        auto(:,4) = auto(:,4) + mHist(tmpz2, autotime, double(ind2(1:end-1)==1).*double(ind2(2:end)==0));

                        subplot('position', [0.875 0.2 0.05 0.6]);
                        bar(0,cnt1/t_Counts);
                        axis([-0.4 0.4 0 1]);
                        set(gca,'xtick',[],'ytick',[]);
                        subplot('position', [0.1 0.1 0.7 0.8]);
                        plot(autotime(1:end-1)*timeunit*1e-9,auto(1:end-1,:)/sum(time)/timeunit*1e9);
                        drawnow;
                    end;
                end;
                n = n+1;
            end

        end
        clf

        autotime = autotime(1:end-1)*timeunit*1e-9;
        auto = auto(1:end-1,:)/sum(time)/timeunit*1e9;

        plot(-autotime(end:-1:2),auto(end:-1:2,2),autotime,auto(:,1))

        res.autotime = autotime;
        res.auto     = auto;
        res.tau      = tau(sum(tcspc,2)>0)*Resolution;
        res.tcspc    = tcspc(sum(tcspc,2)>0,:);
        res.rate     = rate;
        res.time     = time;
    end


    if strcmp(taskflag,'nofcs')

        if nargin<3
            Nsub     = 10;
            timeunit = sync;
            mct      = 3;          % 3 sec max correlation time
            nanz     = 1;
        else
            mct   = param(1);
            Nsub  = param(2);
            if length(param)==2
                timeunit = Resolution;
            else
                timeunit = param(3)*Resolution;
            end
        end
        if nargin == 5
            nanz = weight;
        end;

        Ncasc    = ceil(log2(mct/(timeunit*1e-9*Nsub)+1));
        autotime = 1e-9*timeunit*cumsum(reshape(repmat(2.^[0 0:Ncasc-1],Nsub,1),(Ncasc+1)*Nsub,1));
        auto     = zeros(length(autotime),2);

        timetmp  = ((ones(nanz,1).*1e9/timeunit).^cumsum(ones(nanz,1)))*ones(1,length(autotime));

        rate    = [];
        time    = [];
        tst     = 1;
        cnt     = 0;
        cnt1    = 0;

        max_time  = 0;
        overcount = 0;

        while (tst>0)

            max_time = max_time + timeframe;
            if end_time > 0
                if max_time > end_time
                    max_time = end_time;
                end;
            end;

            [tmpy, tmpx, flag, tst, overcount] = ingopt3Read_t(name, cnt1+1, overcount, max_time, sync);

            if (tst>0)

                cnt1 = cnt1 + tst;

                tmpy  = round((tmpy+tmpx*Resolution)/timeunit);  % in timeunits
                tmpy  = tmpy-tmpy(1);

                cnt  = cnt+length(tmpy);
                dt   = (tmpy(end))*timeunit*1e-9;
                rate = [rate length(tmpy)/dt];
                time = [time dt];

                if tmpy(end)>z32
                    auto = auto +  ingopt32xfcs(tmpy, double(flag<=2), Ncasc, Nsub, nanz);
                else
                    auto = auto + tttr2xfcs(tmpy, double(flag<=2), Ncasc, Nsub, nanz);
                end

                subplot('position', [0.95 0.15 0.025 0.7]);
                bar(0, cnt1/t_Counts);
                axis([-0.4 0.4 0 1]);
                set(gca,'xtick',[],'ytick',[]);
                subplot('position', [0.1 0.1 0.8 0.8]);
                loglog(autotime,reshape(auto,size(autotime,1),nanz).*timetmp'/sum(time),'x');
                drawnow;
            end

        end
        clf

        auto = reshape(auto,size(autotime,1),nanz).*timetmp'/sum(time);
        loglog(autotime,sum(auto,3),'x');

        res.autotime = autotime;
        res.auto     = auto;
        res.rate     = rate;
        res.time     = time;
    end


    if strcmp(taskflag,'nfcs')

        res0    = ingopt3Pro(name, 'tcspc');                %   get tcspc histogram
        tau     = res0.tau;
        tcspc   = res0.tcspc;

        [lt, c] = tcspc_fit(1:length(tau), tcspc, 1, Resolution);            %   fit single exponential

        t1 = 1:length(tau);
        t1 = t1(tau<=20);
        t2 = t1;

        t1 = t1(tcspc(:,1)==max(tcspc(:,1)))+25:t1(end);
        t2 = t2(tcspc(:,2)==max(tcspc(:,2)))+25:t2(end);

        if length(t1)<length(t2)
            t2 = t2(1:length(t1));
        else
            t1 = t1(1:length(t2));
        end;

        [err, ccy1, zzcy1] = ExpFun(lt(1), tau(t1), tcspc(t1,1)-c(1,1));
        [err, ccy2, zzcy2] = ExpFun(lt(2), tau(t2), tcspc(t2,2)-c(2,1));

        tau1 = tau(t1)';
        tau2 = tau(t2)';

        % filter generation

        weight(:,1) = tcspc(tau>=tau1(1) & tau<=tau1(end),1);
        weight(:,2) = tcspc(tau>=tau2(1) & tau<=tau2(end),2);

        para   = [zzcy1(:,1) zzcy2(:,1) ones(size(zzcy1,1),2)];

        for j=1:size(para,2)
            para(:,j) = para(:,j)/sum(para(:,j));     % normalization
        end

        master =          (para(:,[1 3])'*diag(1./weight(:,1))*para(:,[1 3]))\(diag(1./weight(:,1))*para(:,[1 3]))';
        master = [master; (para(:,[2 4])'*diag(1./weight(:,2))*para(:,[2 4]))\(diag(1./weight(:,2))*para(:,[2 4]))']';

        end_time = 0;

        if nargin<3
            Nsub     = 10;
            timeunit = sync;
            mct      = 3;          % 3 sec max correlation time
            nanz     = 1;
        else
            mct   = param(1);
            Nsub  = param(2);
            if length(param)==2
                timeunit = Resolution;
            else
                timeunit = param(3)*Resolution;
            end
        end
        if nargin == 4
            nanz = weight;
        end;

        Ncasc    = ceil(log2(mct/(timeunit*1e-9*Nsub)+1));
        autotime = 1e-9*timeunit*cumsum(reshape(repmat(2.^[0 0:Ncasc-1],Nsub,1),(Ncasc+1)*Nsub,1));
        auto     = zeros(nanz*length(autotime),1);

        timetmp  = ((ones(nanz,1).*1e9/timeunit).^cumsum(ones(nanz,1)))*ones(1,length(autotime));

        tau      = (0:(NChannels-1));
        tcspc    = zeros(NChannels,2);

        rate    = [];
        time    = [];
        tst1    = 1;
        cnt     = 0;
        cnt1    = 0;
        scnt    = 0;

        overcount = 0;

        max_time  = 0;
        intensity = 0;

        while (tst1>0)

            max_time = max_time + timeframe;
            if end_time > 0
                if max_time > end_time
                    max_time = end_time;
                end;
            end;

            [tmpy, tmpx, flag, tst, overcount] = ingopt3Read_t(name, cnt1+1, overcount, max_time, sync);

            if (tst1>0)

                cnt1 = cnt1 + tst;
                scnt = scnt + 1;

                tmp1  = (flag==1)&&(tmpx>=tau1(1))&&(tmpx<=tau1(end));  % select valid tcspc-times
                tmp2  = (flag==2)&&(tmpx>=tau2(1))&&(tmpx<=tau2(end));  % select valid tcspc-times

                tmpy(~(tmp1||tmp2)) = [];
                tmpx(~(tmp1||tmp2)) = [];
                flag(~(tmp1||tmp2)) = [];

                tmpy  = round((tmpy+tmpx*Resolution)/timeunit);  % in timeunits
                tmpy  = tmpy-tmpy(1);

                tmpx(tmp1) = tmpx(tmp1)-tau1(1) + 1;
                tmpx(tmp2) = tmpx(tmp2)-tau2(1) + 1;

                cnt  = cnt+length(tmpy);
                dt   = (tmpy(end))*timeunit*1e-9;
                rate = [rate length(tmpy)/dt];
                time = [time dt];

                imp  = master(tmpx,1).*(flag==1)+master(tmpx,3).*(flag==2);

                intensity(scnt) = sum(imp);

                if tmpy(end)>z32
                    auto = auto +  ingopt32xfcs(tmpy, imp, Ncasc, Nsub, nanz);
                else
                    auto = auto + tttr2xfcs(tmpy, imp, Ncasc, Nsub, nanz);
                end

                subplot('position', [0.875 0.2 0.05 0.6]);
                bar(0, cnt1/t_Counts);
                axis([-0.4 0.4 0 1]);
                set(gca,'xtick',[],'ytick',[]);
                subplot('position', [0.1 0.1 0.7 0.8]);
                semilogx(autotime,reshape(auto,size(autotime,1),nanz).*timetmp'/sum(time));
                drawnow;
            end

        end
        clf

        auto = reshape(auto,size(autotime,1),nanz).*timetmp'/sum(time);
        loglog(autotime,sum(auto,3));

        res.autotime  = autotime;
        res.auto      = auto;
        res.intensity = intensity;
        res.tau       = Resolution*tau(sum(tcspc,2)>0);
        res.tcspc     = tcspc(sum(tcspc,2)>0,:);
        res.rate      = rate;
        res.time      = time;

    end





    if strcmp(taskflag,'hofcs')  % higher-order correlation
        if nargin<3
            Nsub  = 10;
            timeunit = 1e4;
            Ncasc = ceil(log2(3/(timeunit*1e-9*Nsub)+1));    % 3 sec max correlation time
            Nord = 4;
        else
            Nsub  = 10;
            timeunit = param(1);
            Ncasc = ceil(log2(3/(timeunit*1e-9*Nsub)+1));    % 3 sec max correlation time
            Nord = param(2);
        end

        autotime = 1e-9*timeunit*cumsum(reshape(repmat(2.^[0 0:Ncasc-1],Nsub,1),(Ncasc+1)*Nsub,1));
        auto  = zeros(length(autotime),Nord,Nord,2);
        intensity = zeros(Nord,2);

        tau   = (0:(NChannels-1));
        tcspc = zeros(NChannels,2);

        rate    = [];
        time    = [];
        tst     = 1;
        cnt     = 0;
        cnt1    = 0;

        max_time  = 0;
        overcount = 0;

        while (tst>0)

            max_time = max_time + timeframe;

            [tmpy, tmpx, flag, tst, overcount] = ingopt3Read_t(name, cnt1+1, overcount, max_time, sync);

            if (tst>0)

                cnt1 = cnt1 + tst;

                tmpy = round((tmpy+tmpx*Resolution)/timeunit);  % in timeunits
                tmpy = tmpy-tmpy(1);

                cnt = cnt+length(tmpy);
                dt = (tmpy(end))*timeunit*1e-9;
                rate = [rate length(tmpy)/dt];
                time = [time dt];

                tcspc(:,1) = tcspc(:,1) + mHist(tmpx(flag==1), tau);
                tcspc(:,2) = tcspc(:,2) + mHist(tmpx(flag==2), tau);

                tmpy = ceil(tmpy);
                ind1 = double(flag==1);
                ind2 = double(flag==2);

                t = 1:length(tmpy);
                dtmp = diff(tmpy);
                t1 = t([1; dtmp]>0 & [dtmp; 1]==0);
                t2 = t([1; dtmp]==0 & [dtmp; 1]>0);
                for j=1:length(t1)
                    ind1(t1(j)) = sum(ind1(t1(j):t2(j)));
                    ind2(t1(j)) = sum(ind2(t1(j):t2(j)));
                end
                ind = [1; dtmp]==0;
                tmpy(ind) = [];
                ind1(ind) = [];
                ind2(ind) = [];

                for j=1:Nord
                    intensity(j,1) = intensity(j,1) + sum(ind1.^j);
                    intensity(j,2) = intensity(j,2) + sum(ind2.^j);
                    for k=1:Nord
                        if tmpy(end)>z32
                            auto(:,j,k,1) = auto(:,j,k,1) + ingopt32fcs(tmpy, ind1.^j, ind2.^k, Ncasc, Nsub, smooth);
                            auto(:,j,k,2) = auto(:,j,k,2) + ingopt32fcs(tmpy, ind2.^j, ind1.^k, Ncasc, Nsub, smooth);
                        else
                            auto(:,j,k,1) = auto(:,j,k,1) + tttr2fcs(tmpy, ind1.^j, ind2.^k, Ncasc, Nsub, smooth);
                            auto(:,j,k,2) = auto(:,j,k,2) + tttr2fcs(tmpy, ind2.^j, ind1.^k, Ncasc, Nsub, smooth);
                        end
                    end
                end

                subplot('position', [0.875 0.2 0.05 0.6]);
                bar(0, cnt1/t_Counts);
                axis([-0.4 0.4 0 1]);
                set(gca,'xtick',[],'ytick',[]);
                subplot('position', [0.1 0.1 0.7 0.8]);
                bld = auto(:,1,1,1) + auto(:,1,1,2);
                for j=2:Nord
                    bld = [bld auto(:,j,j,1) + auto(:,j,j,2)];
                end
                semilogx(autotime,bld/sum(time)/timeunit*1e3);
                drawnow;
            end
        end
        clf

        bld = auto(:,1,1,1) + auto(:,1,1,2);
        for j=2:Nord
            bld = [bld auto(:,j,j,1) + auto(:,j,j,2)];
        end
        semilogx(autotime, bld/sum(time)/timeunit*1e3);

        res.autotime  = autotime;
        res.auto      = auto/sum(time)/timeunit*1e3;
        res.intensity = intensity/sum(time);
        res.tau   = tau(sum(tcspc,2)>0)*Resolution;
        res.tcspc = tcspc(sum(tcspc,2)>0,:);
        res.rate = rate;
        res.time = time;
    end


    if strcmp(taskflag,'mlcs')
        if nargin<3
            disp('No master patterns given!');
            return
        end

        Nsub  = 10;
        Ncasc = ceil(log2(3/(sync*1e-9*Nsub)+1));    % 3 sec max correlation time
        autotime = 1e-9*sync*cumsum(reshape(repmat(2.^[0 0:Ncasc-1],Nsub,1),(Ncasc+1)*Nsub,1));
        timeunit = sync;

        tcspc = zeros(size(param,1),2);

        rate    = [];
        time    = [];
        tst     = 1;
        cnt     = 0;
        cnt1    = 0;

        max_time  = 0;
        overcount = 0;

        tau = param(:,1);
        param = param(:,2:end);

        % ML filter generation
        % para = [dye1 channel#1; dye1 channel#2; dye2 channel#1; dye2 channel#2];

        for j=1:size(param,2)
            param(:,j) = param(:,j)/sqrt(sum(param(:,j).^2));  % normalization
        end

        master = [param(:,1:2:end) param(:,2:2:end)];
        len = size(master,2)/2;

        auto = zeros(length(autotime),len,len);

        while (tst>0)

            max_time = max_time + timeframe;

            [tmpy, tmpx, flag, tst, overcount] = ingopt3Read_t(name, cnt1+1, overcount, max_time, sync);

            if (tst>0)

                cnt1 = cnt1 + tst;
                tmpx = tmpx*Resolution;             % in ns
                tmpy = round((tmpy+tmpx)/timeunit);  % in timeunits

                tst = (tmpx>=tau(1)&tmpx<=tau(end));
                tmpy = tmpy(tst);
                tmpx = tmpx(tst);
                flag = flag(tst);

                tcspc(:,1) = tcspc(:,1) + mHist(tmpx(flag==1), tau);
                tcspc(:,2) = tcspc(:,2) + mHist(tmpx(flag==2), tau);

                tmpy = tmpy-tmpy(1);
                tmpx = 1+ round((tmpx-tau(1))/Resolution);

                cnt  = cnt+length(tmpy);
                dt   = (tmpy(end))*timeunit*1e-9;
                rate = [rate length(tmpy)/dt];
                time = [time dt];

                for j=1:len
                    for k=1:len
                        if tmpy(end)>z32
                            auto(:,j,k) = auto(:,j,k) + ingopt32fcs(tmpy, master(tmpx,j).*double(flag==1), master(tmpx,end/2+k).*double(flag==2), Ncasc, Nsub, smooth);
                            auto(:,j,k) = auto(:,j,k) + ingopt32fcs(tmpy, master(tmpx,end/2+k).*double(flag==2), master(tmpx,j).*double(flag==1), Ncasc, Nsub, smooth);
                        else
                            auto(:,j,k) = auto(:,j,k) + tttr2fcs(tmpy, master(tmpx,j).*double(flag==1), master(tmpx,end/2+k).*double(flag==2), Ncasc, Nsub, smooth);
                            auto(:,j,k) = auto(:,j,k) + tttr2fcs(tmpy, master(tmpx,end/2+k).*double(flag==2), master(tmpx,j).*double(flag==1), Ncasc, Nsub, smooth);
                        end
                    end
                end

                subplot('position', [0.925 0.2 0.025 0.6]);
                bar(0, cnt1/t_Counts);
                axis([-0.4 0.4 0 1]);
                set(gca,'xtick',[],'ytick',[]);
                subplot('position', [0.1 0.1 0.7 0.8]);
                semilogx(autotime, auto(:,:,1)/sum(time)/timeunit*1e3);
                drawnow;
            end
        end
        clf
        semilogx(autotime,auto(:,:,1)/sum(time)/timeunit*1e3);

        auto = auto/sum(time)/timeunit*1e3;

        res.autotime = autotime;
        res.auto0    = auto;

        if 0
            for j=1:len
                for jj=1:len
                    for k=1:len
                        for kk=1:len
                            mm((j-1)*len+k,(jj-1)*len+kk) = (master(:,j)'*master(:,end/2+jj))*(master(:,k)'*master(:,end/2+kk)) + ...
                                (master(:,end/2+j)'*master(:,jj))*(master(:,end/2+k)'*master(:,kk));
                        end
                    end
                end
            end
            res.mm = mm;
            auto = permute(auto,[2 3 1]);
            for s=1:size(auto,3)
                auto(:,:,s) = reshape(mm\reshape(squeeze(auto(:,:,s)),len^2,1),len,len);
            end
            auto = permute(auto,[3 1 2]);
            res.auto = auto;
        end

        res.tau = tau;
        res.tcspc = tcspc;
        res.rate = rate;
        res.master = master;
        res.time = time;
    end



    if strcmp(taskflag,'fimda')
        Ncasc = 10;
        Nsub = 10;
        nums = 10000;
        timeunit = sync;

        autotime = 1e-9*timeunit*cumsum(reshape(repmat(2.^(0:Ncasc-1),Nsub,1),Ncasc*Nsub,1));
        fimda    = zeros(nums+1,length(autotime));

        tau   = (0:(NChannels-1));
        tcspc = zeros(NChannels,2);

        rate    = [];
        time    = [];
        tst     = 1;
        cnt     = 0;
        cnt1    = 0;

        max_time  = 0;
        overcount = 0;

        while (tst>0)

            max_time = max_time + timeframe;

            [tmpy, tmpx, flag, tst, overcount] = ingopt3Read_t(name, cnt1+1, overcount, max_time, sync);

            if (tst>0)

                cnt1 = cnt1 + tst;

                tmpy = round((tmpy+tmpx*Resolution)/timeunit);  % in timeunits
                cnt = cnt+length(tmpy);
                dt = (tmpy(end) - tmpy(1))*timeunit;
                rate = [rate length(tmpy)/dt];
                time = [time dt];

                tcspc(:,1) = tcspc(:,1) + mHist(tmpx(flag==1), tau);
                tcspc(:,2) = tcspc(:,2) + mHist(tmpx(flag==2), tau);

                fimda = fimda + tttr2fimda(tmpy, ones(size(tmpy)), Ncasc, Nsub, nums);
                subplot('position', [0.875 0.2 0.05 0.6]);
                bar(0, cnt1/t_Counts);
                axis([-0.4 0.4 0 1]);
                set(gca,'xtick',[],'ytick',[]);
                subplot('position', [0.1 0.1 0.7 0.8]);
                loglog(0:nums,fimda);
                drawnow;
            end

        end
        clf

        loglog(0:nums,fimda);

        res.nums     = 0:nums;
        res.fimda    = fimda;
        res.autotime = autotime;
        res.tau      = Reaolution*tau(sum(tcspc,2)>0);
        res.tcspc    = tcspc(sum(tcspc,2)>0,:);
        res.rate     = rate;
        res.time     = time;

    end

    head.NCounts  = cnt;
end;
warning on MATLAB:DeprecatedLogicalAPI