function [res, head] = pt3Pro(name1, name2, flag, param, weight);

warning off MATLAB:DeprecatedLogicalAPI

timeframe = 6e10; % time-frame to read photons from file in ns
smooth = 1; % smoothed autocorrelation
close all   

head  = pt3Read_t(name1);
head2 = pt3Read_t(name2);
CntRate0 = (head.CntRate0+head2.CntRate0)/2;


sync = 1e9/CntRate0; 

if strcmp(flag,'tcspc')
    z32 = 2^32-1;

    NChannels = floor(sync/head.Resolution);
    tau   = (0:(NChannels-1))*head.Resolution;
    tcspc = zeros(NChannels,2);
    
    head = pt3Read_t(name1);   
    t_Counts = head.NCounts;
    head = pt3Read_t(name2);        
    t_Counts = t_Counts + head.NCounts;
    
    rate    = [];
    time    = [];
    tst1    = 1; 
    tst2    = 1; 
    cnt     = 0;
    cnt1    = 0; 
    cnt2    = 0; 
 
    max_time  = 0;

    overcount1 = 0;
    overcount2 = 0;
    
    while ((tst1>0)&(tst2>0))

        max_time = max_time + timeframe;

        [head, tmpy1, tmpx1, z1, tst1, overcount1] = pt3Read_t(name1, cnt1+1, overcount1, max_time, CntRate0);
        [head, tmpy2, tmpx2, z2, tst2, overcount2] = pt3Read_t(name2, cnt2+1, overcount2, max_time, CntRate0);

        cnt1 = cnt1 + tst1;
        cnt2 = cnt2 + tst2;

        if ((tst1>0)&(tst2>0))

            tcspc(:,1) = tcspc(:,1) + mHist(tmpx1, tau);
            tcspc(:,2) = tcspc(:,2) + mHist(tmpx2, tau);

            subplot('position', [0.875 0.2 0.05 0.6]);
            bar(0,(cnt1+cnt2)/t_Counts);
            axis([-0.4 0.4 0 1]);
            set(gca,'xtick',[],'ytick',[]);
            subplot('position', [0.1 0.1 0.7 0.8]);
            semilogy(tau,tcspc);
            drawnow;
        end

    end
    clf

    semilogy(tau,tcspc);
    res = struct('tau', tau);
    res = setfield(res, 'tcspc', tcspc);
end


if strcmp(flag,'bin')

    timeunit = sync;
    if nargin<4
        binwidth = 1e5/timeunit; % in timeunits
    else
        binwidth = param(1)/timeunit; % in timeunits
    end
    
    head = pt3Read_t(name1);   
    t_Counts = head.NCounts;
    head = pt3Read_t(name2);        
    t_Counts = t_Counts + head.NCounts;
    
    rate1   = [];
    time1   = [];
    bin1    = [];    
    rate2   = [];
    time2   = [];
    bin2    = [];
    tst1    = 1; 
    tst2    = 1; 
    cnt     = 0;
    cnt1    = 0; 
    cnt2    = 0; 
 
    max_time  = 0;

    overcount1 = 0;
    overcount2 = 0;
    
    while ((tst1>0)&(tst2>0))

        max_time = max_time + timeframe;

        [head, tmpy1, tmpx1, z1, tst1, overcount1] = pt3Read_t(name1, cnt1+1, overcount1, max_time, CntRate0);
        [head, tmpy2, tmpx2, z2, tst2, overcount2] = pt3Read_t(name2, cnt2+1, overcount2, max_time, CntRate0);

        cnt1 = cnt1 + tst1;
        cnt2 = cnt2 + tst2;

        if ((tst1>0)&(tst2>0))

            tmpy1 = round([tmpy1+tmpx1]/timeunit);  % in timeunits
            tmpy2 = round([tmpy2+tmpx2]/timeunit);  % in timeunits

            tmpy1 = tmpy1-tmpy1(1);
            tmpy2 = tmpy2-tmpy2(1);

            dt1   = tmpy1(end)*timeunit*1e-9;
            rate1 = [rate1 length(tmpy1)/dt1];
            time1 = [time1 dt1];

            dt2   = tmpy2(end)*timeunit*1e-9;
            rate2 = [rate2 length(tmpy2)/dt2];
            time2 = [time2 dt2];

            if tmpy1(end)>tmpy2(end)
                tmpy2(end+1)=tmpy1(end); 
                bin1 = [bin1 tttr2bin(tmpy1, binwidth)'];
                bin2 = [bin2 tttr2bin(tmpy2, binwidth)'];
                bin2(end) = bin2(end)-1;
            else
                tmpy1(end+1)=tmpy2(end); 
                bin1 = [bin1 tttr2bin(tmpy1, binwidth)'];
                bin2 = [bin2 tttr2bin(tmpy2, binwidth)'];
                bin1(end) = bin1(end)-1;                
            end

            subplot('position', [0.875 0.2 0.05 0.6]);
            bar(0,(cnt1+cnt2)/t_Counts);
            axis([-0.4 0.4 0 1]);
            set(gca,'xtick',[],'ytick',[]);
            subplot('position', [0.1 0.1 0.7 0.8]);
            plot(1:1e4,bin1(end-1e4+1:end));
            axis tight
            drawnow;
        end

    end

    res = struct('bin1', bin1);
	res = setfield(res, 'bin2', bin2);
    res = setfield(res, 'binwidth', binwidth*timeunit*1e-9);
    res = setfield(res, 'rate1', rate1);
    res = setfield(res, 'time1', time1);
    res = setfield(res, 'rate2', rate2);
    res = setfield(res, 'time2', time2);
end


if strcmp(flag,'fcs')
    end_time = 0;
    if nargin<4
        Nsub  = 10;
        Ncasc = ceil(log2(3/(sync*1e-9*Nsub)+1));    % 3 sec max correlation time
        autotime = 1e-9*sync*cumsum(reshape(repmat(2.^[0 0:Ncasc-1],Nsub,1),(Ncasc+1)*Nsub,1));
        auto  = zeros(length(autotime),3);
        timeunit = sync;
    else
        Ncasc = param(1);
        Nsub = param(2);
        if length(param)==2
            timeunit = head.Resolution;
        else
            timeunit = param(3)*head.Resolution;
            if length(param)>3
                end_time = param(4)*1e9;
            end;
        end
        autotime = 1e-9*timeunit*cumsum(reshape(repmat(2.^[0 0:Ncasc-1],Nsub,1),(Ncasc+1)*Nsub,1));
        auto     = zeros(length(autotime),3);
    end
    
    z32 = 2^32-1;

    NChannels = floor(sync/head.Resolution);
    tau   = (0:(NChannels-1))*head.Resolution;
    tcspc = zeros(NChannels,2);
    
    head = pt3Read_t(name1);   
    t_Counts = head.NCounts;
    head = pt3Read_t(name2);        
    t_Counts = t_Counts + head.NCounts;
    
    rate    = [];
    time    = [];
    tst1    = 1; 
    tst2    = 1; 
    cnt     = 0;
    cnt1    = 0; 
    cnt2    = 0; 
 
    max_time  = 0;

    overcount1 = 0;
    overcount2 = 0;
    
    while ((tst1>0)&(tst2>0))

        max_time = max_time + timeframe;
        if end_time > 0 
            if max_time > end_time
                max_time = end_time;
            end;
        end;

        [head, tmpy1, tmpx1, z1, tst1, overcount1] = pt3Read_t(name1, cnt1+1, overcount1, max_time, CntRate0);
        [head, tmpy2, tmpx2, z2, tst2, overcount2] = pt3Read_t(name2, cnt2+1, overcount2, max_time, CntRate0);

        cnt1 = cnt1 + tst1;
        cnt2 = cnt2 + tst2;

        if ((tst1>0)&(tst2>0))

            tmpy = round([tmpy1+tmpx1; tmpy2+tmpx2]/timeunit);  % in timeunits
            tmpx = [tmpx1; tmpx2];                              % in ns

            ind  = [zeros(length(tmpy),1)];
            ind1 = [zeros(length(tmpy1),1); ones(length(tmpy2),1)];

            [tmpy, ord] = sort(tmpy);
            tmpx = tmpx(ord);
            ind1 = ind1(ord);
            tmpy = tmpy-tmpy(1);

            cnt = cnt+length(tmpy);
            dt = (tmpy(end))*timeunit*1e-9;  % in s
            rate = [rate length(tmpy)/dt];
            time = [time dt];

            tcspc(:,1) = tcspc(:,1) + mHist(tmpx1, tau);
            tcspc(:,2) = tcspc(:,2) + mHist(tmpx2, tau);

            if tmpy(end)>z32
                auto(:,1) = auto(:,1) + pt32fcs(tmpy, double(ind1==0), double(ind1==0), Ncasc, Nsub, smooth);
                auto(:,2) = auto(:,2) + pt32fcs(tmpy, double(ind1==1), double(ind1==1), Ncasc, Nsub, smooth);
                auto(:,3) = auto(:,3) + pt32fcs(tmpy, double(ind==0), double(ind==0), Ncasc, Nsub, smooth);
            else
                auto(:,1) = auto(:,1) + tttr2fcs(tmpy, double(ind1==0), double(ind1==0), Ncasc, Nsub, smooth);
                auto(:,2) = auto(:,2) + tttr2fcs(tmpy, double(ind1==1), double(ind1==1), Ncasc, Nsub, smooth);
                auto(:,3) = auto(:,3) + tttr2fcs(tmpy, double(ind==0), double(ind==0), Ncasc, Nsub, smooth);
            end

            subplot('position', [0.875 0.2 0.05 0.6]);
            bar(0,(cnt1+cnt2)/t_Counts);
            axis([-0.4 0.4 0 1]);
            set(gca,'xtick',[],'ytick',[]);
            subplot('position', [0.1 0.1 0.7 0.8]);
            semilogx(autotime,auto/sum(time)/timeunit*1e9);
            drawnow;
        end

    end
    clf

    semilogx(autotime,auto/sum(time)/timeunit*1e9);
    res = struct('autotime', autotime);
    res = setfield(res, 'auto', auto/sum(time)/timeunit*1e9);
    tmp = 0:head.NChannels;
    res = setfield(res, 'tau', tau);
    res = setfield(res, 'tcspc', tcspc);
    res = setfield(res, 'rate', rate);
    res = setfield(res, 'time', time);
    
end


if strcmp(flag,'nfcs')
    
    [res0, head]  = pt3Pro(name1, name2, 'tcspc');                %   get tcspc histogram
    tau           = res0.tau;
    tcspc         = res0.tcspc;

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

    lt1 = simplex('ExpFun',[1 5],[],[],[],[],tau(t1),tcspc(t1,1));
    lt2 = simplex('ExpFun',[1 5],[],[],[],[],tau(t2),tcspc(t2,2));    
    [err, ccy1] = ExpFun(lt1, tau(t1), tcspc(t1,1));
    [err, ccy2] = ExpFun(lt2, tau(t2), tcspc(t2,2));
    ccy1(1) = 0; ccy2(1) = 0; 
    zzcy1 = ExpFun(lt1, ccy1, tau(t1));
    zzcy2 = ExpFun(lt2, ccy2, tau(t2));
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

    z32      = 2^32-1;
    end_time = 0;
    
    timeunit = sync;
    Nsub     = 10;
    Ncasc    = ceil(log2(3/(sync*1e-9*Nsub)+1));    % 3 sec max correlation time
    flag     = 5;
    
    if nargin==4
        flag = param(1);
    end
    
    autotime = 1e-9*timeunit*cumsum(reshape(repmat(2.^[0 0:Ncasc-1],Nsub,1),(Ncasc+1)*Nsub,1));
    auto     = zeros(flag*length(autotime),1);
    
    timetmp  = ((ones(flag,1).*1e9/timeunit).^cumsum(ones(flag,1)))*ones(1,length(autotime));    

    NChannels = floor(sync/head.Resolution);
    tau   = (0:(NChannels-1))*head.Resolution;
    tcspc = zeros(NChannels,2);
    
    head = pt3Read_t(name1);   
    t_Counts = head.NCounts;
    head = pt3Read_t(name2);        
    t_Counts = t_Counts + head.NCounts;
        
    rate    = [];
    time    = [];
    tst1    = 1; 
    tst2    = 1; 
    cnt     = 0;
    cnt1    = 0; 
    cnt2    = 0; 
    scnt = 0;
 
    max_time  = 0;

    overcount1 = 0;
    overcount2 = 0;          

    intensity = 0;
    while ((tst1>0)&(tst2>0))

        max_time = max_time + timeframe;
        if end_time > 0 
            if max_time > end_time
                max_time = end_time;
            end;
        end;

        [head, tmpy1, tmpx1, z1, tst1, overcount1] = pt3Read_t(name1, cnt1+1, overcount1, max_time, CntRate0);
        [head, tmpy2, tmpx2, z2, tst2, overcount2] = pt3Read_t(name2, cnt2+1, overcount2, max_time, CntRate0);

        cnt1 = cnt1 + tst1;
        cnt2 = cnt2 + tst2;

        if ((tst1>0)&(tst2>0))

            scnt = scnt + 1;
            
            ind1  = zeros(length(tmpy1),1); 
            ind2  =  ones(length(tmpy2),1);
            
            tmp1  = (tmpx1>=tau1(1)) & (tmpx1<=tau1(end));         % select valid tcspc-times
            tmp2  = (tmpx2>=tau2(1)) & (tmpx2<=tau2(end));         % select valid tcspc-times

            tmpy1 = tmpy1(tmp1);
            tmpy2 = tmpy2(tmp2);
            tmpx1 = tmpx1(tmp1);
            tmpx2 = tmpx2(tmp2);
            ind1  = ind1(tmp1);
            ind2  = ind2(tmp2);


            ind   = [ind1; ind2];
            tmpy  = round([tmpy1+tmpx1; tmpy2+tmpx2]/timeunit);  % in timeunits

            tmpx1 = round((tmpx1-tau1(1))/head.Resolution)+1;
            tmpx2 = round((tmpx2-tau2(1))/head.Resolution)+1;
            tmpx  = [tmpx1; tmpx2];                              % in ns

            [tmpy, ord] = sort(tmpy);                            % sort arrival times in asc. order
            tmpx  = tmpx(ord);
            ind   = ind(ord);
            tmpy  = tmpy-tmpy(1);
            
            cnt  = cnt+length(tmpy);
            dt   = (tmpy(end))*timeunit*1e-9;
            rate = [rate length(tmpy)/dt];
            time = [time dt];

            imp  = master(tmpx,1).*(ind==0)+master(tmpx,3).*(ind==1);

            intensity(scnt) = sum(imp);
            
            if tmpy(end)>z32
                auto(:,scnt) = pt32nfcs(tmpy, imp, Ncasc, Nsub, flag);
            else
                auto(:,scnt) = tttr2nfcs(tmpy, imp, Ncasc, Nsub, flag);
            end

            subplot('position', [0.875 0.2 0.05 0.6]);
            bar(0,(cnt1+cnt2)/t_Counts);
            axis([-0.4 0.4 0 1]);
            set(gca,'xtick',[],'ytick',[]);
            subplot('position', [0.1 0.1 0.7 0.8]);
            loglog(autotime,reshape(sum(auto,2),size(autotime,1),flag).*timetmp'/sum(time));
            drawnow;
        end

    end
    clf

            
    auto = reshape(auto.*repmat(reshape(timetmp',size(auto,1),1),[1 size(auto,2)]),[size(autotime,1),flag,size(auto,2)])/sum(time);
    loglog(autotime,sum(auto,3));
    res = struct('autotime', autotime);
    res = setfield(res, 'auto', auto);
    res = setfield(res, 'intensity', intensity);
    tmp = 0:head.NChannels;
    res = setfield(res, 'tau', tau);
    res = setfield(res, 'tcspc', tcspc);
    res = setfield(res, 'rate', rate);
    res = setfield(res, 'time', time);
    
end


if strcmp(flag,'xfcs')
    if nargin<4
        Nsub  = 10;
        Ncasc = ceil(log2(3/(sync*1e-9*Nsub)+1));    % 3 sec max correlation time
        autotime = 1e-9*sync*cumsum(reshape(repmat(2.^[0 0:Ncasc-1],Nsub,1),(Ncasc+1)*Nsub,1));
        auto  = zeros(length(autotime),2);
        timeunit = sync;
    else
        Ncasc = param(1);
        Nsub = param(2);
        if length(param)==2
            timeunit = head.Resolution;
        else
            timeunit = param(3)*head.Resolution;
        end
        autotime = 1e-9*timeunit*cumsum(reshape(repmat(2.^[0 0:Ncasc-1],Nsub,1),(Ncasc+1)*Nsub,1));
        auto  = zeros(length(autotime),2);
    end
    
    z32 = 2^32-1;

    NChannels = floor(sync/head.Resolution);
    tau   = (0:(NChannels-1))*head.Resolution;
    tcspc = zeros(NChannels,2);
    
    head = pt3Read_t(name1);   
    t_Counts = head.NCounts;
    head = pt3Read_t(name2);        
    t_Counts = t_Counts + head.NCounts;
    
    rate    = [];
    time    = [];
    tst1    = 1; 
    tst2    = 1; 
    cnt     = 0;
    cnt1    = 0; 
    cnt2    = 0; 
 
    max_time  = 0;

    overcount1 = 0;
    overcount2 = 0;
    
    while ((tst1>0)&(tst2>0))

        max_time = max_time + timeframe;

        [head, tmpy1, tmpx1, z1, tst1, overcount1] = pt3Read_t(name1, cnt1+1, overcount1, max_time, CntRate0);
        [head, tmpy2, tmpx2, z2, tst2, overcount2] = pt3Read_t(name2, cnt2+1, overcount2, max_time, CntRate0);

        cnt1 = cnt1 + tst1;
        cnt2 = cnt2 + tst2;

        if ((tst1>0)&(tst2>0))

            tmpy = round([tmpy1+tmpx1; tmpy2+tmpx2]/timeunit);  % in timeunits
            tmpx = [tmpx1; tmpx2];                                  % in ns
            ind = [zeros(length(tmpy1),1); ones(length(tmpy2),1)];

            [tmpy, ord] = sort(tmpy);
            tmpx = tmpx(ord);
            ind = ind(ord);
            tmpy = tmpy-tmpy(1);

            cnt = cnt+length(tmpy);
            dt = (tmpy(end))*timeunit*1e-9;
            rate = [rate length(tmpy)/dt];
            time = [time dt];

            tcspc(:,1) = tcspc(:,1) + mHist(tmpx1, tau);
            tcspc(:,2) = tcspc(:,2) + mHist(tmpx2, tau);

            if tmpy(end)>z32
                auto(:,1) = auto(:,1) + pt32fcs(tmpy, double(ind==0), double(ind==1), Ncasc, Nsub, smooth);
                auto(:,2) = auto(:,2) + pt32fcs(tmpy, double(ind==1), double(ind==0), Ncasc, Nsub, smooth);
            else
                auto(:,1) = auto(:,1) + tttr2fcs(tmpy, double(ind==0), double(ind==1), Ncasc, Nsub, smooth);
                auto(:,2) = auto(:,2) + tttr2fcs(tmpy, double(ind==1), double(ind==0), Ncasc, Nsub, smooth);
            end

            subplot('position', [0.875 0.2 0.05 0.6]);
            bar(0,(cnt1+cnt2)/t_Counts);
            axis([-0.4 0.4 0 1]);
            set(gca,'xtick',[],'ytick',[]);
            subplot('position', [0.1 0.1 0.7 0.8]);
            semilogx(autotime,auto/sum(time)/timeunit*1e9);
            drawnow;
        end

    end
    clf

    semilogx(autotime,auto/sum(time)/timeunit*1e9);
    res = struct('autotime', autotime);
    res = setfield(res, 'auto', auto/sum(time)/timeunit*1e9);
    tmp = 0:head.NChannels;
    res = setfield(res, 'tau', tau);
    res = setfield(res, 'tcspc', tcspc);
    res = setfield(res, 'rate', rate);
    res = setfield(res, 'time', time);
    
end


if strcmp(flag,'xfcsm')
    if nargin<4
        Nsub  = 10;
        Ncasc = ceil(log2(3/(sync*1e-9*Nsub)+1));    % 3 sec max correlation time
        autotime = 1e-9*sync*cumsum(reshape(repmat(2.^[0 0:Ncasc-1],Nsub,1),(Ncasc+1)*Nsub,1));
        auto  = zeros(length(autotime),2);
        timeunit = sync;
    else
        Ncasc = param(1);
        Nsub = param(2);
        if length(param)==2
            timeunit = head.Resolution;
        else
            timeunit = param(3)*head.Resolution;
        end
        autotime = 1e-9*timeunit*cumsum(reshape(repmat(2.^[0 0:Ncasc-1],Nsub,1),(Ncasc+1)*Nsub,1));
        auto  = zeros(length(autotime),2);
    end
    
    z32 = 2^32-1;

    NChannels = floor(sync/head.Resolution);
    tau   = (0:(NChannels-1))*head.Resolution;
    tcspc = zeros(NChannels,2);
    
    head = pt3Read_t(name1);   
    t_Counts = head.NCounts;
    head = pt3Read_t(name2);        
    t_Counts = t_Counts + head.NCounts;
    
    rate    = [];
    time    = [];
    tst1    = 1; 
    tst2    = 1; 
    cnt     = 0;
    cnt1    = 0; 
    cnt2    = 0;
    scnt = 0;
 
    max_time  = 0;

    overcount1 = 0;
    overcount2 = 0;
    
    while ((tst1>0)&(tst2>0))

        max_time = max_time + timeframe;

        [head, tmpy1, tmpx1, z1, tst1, overcount1] = pt3Read_t(name1, cnt1+1, overcount1, max_time, CntRate0);
        [head, tmpy2, tmpx2, z2, tst2, overcount2] = pt3Read_t(name2, cnt2+1, overcount2, max_time, CntRate0);

        cnt1 = cnt1 + tst1;
        cnt2 = cnt2 + tst2;

        if ((tst1>0)&(tst2>0))
            
            scnt = scnt + 1;

            tmpy = round([tmpy1+tmpx1; tmpy2+tmpx2]/timeunit);  % in timeunits

            tmpx = [tmpx1; tmpx2];                                  % in ns
            ind = [zeros(length(tmpy1),1); ones(length(tmpy2),1)];

            [tmpy, ord] = sort(tmpy);
            tmpx = tmpx(ord);
            ind = ind(ord);
            tmpy = tmpy-tmpy(1);

            cnt = cnt+length(tmpy);
            dt = (tmpy(end))*timeunit*1e-9;
            rate = [rate length(tmpy)/dt];
            time = [time dt];

            tcspc(:,1,scnt) = mHist(tmpx1, tau);
            tcspc(:,2,scnt) = mHist(tmpx2, tau);

            if tmpy(end)>z32
                auto(:,1,scnt) = pt32fcs(tmpy, double(ind==0), double(ind==1), Ncasc, Nsub, smooth);
                auto(:,2,scnt) = pt32fcs(tmpy, double(ind==1), double(ind==0), Ncasc, Nsub, smooth);
            else
                auto(:,1,scnt) = tttr2fcs(tmpy, double(ind==0), double(ind==1), Ncasc, Nsub, smooth);
                auto(:,2,scnt) = tttr2fcs(tmpy, double(ind==1), double(ind==0), Ncasc, Nsub, smooth);
            end

            subplot('position', [0.875 0.2 0.05 0.6]);
            bar(0,(cnt1+cnt2)/t_Counts);
            axis([-0.4 0.4 0 1]);
            set(gca,'xtick',[],'ytick',[]);
            subplot('position', [0.1 0.1 0.7 0.8]);
            semilogx(autotime,squeeze(sum(auto,3))/sum(time)/timeunit*1e3);
            drawnow;
        end

    end
    clf

    semilogx(autotime,squeeze(sum(auto,3))/sum(time)/timeunit*1e3);
    res = struct('autotime', autotime);
    res = setfield(res, 'auto', auto/sum(time)/timeunit*1e3);
    tmp = 0:head.NChannels;
    res = setfield(res, 'tau', tau);
    res = setfield(res, 'tcspc', tcspc);
    res = setfield(res, 'rate', rate);
    res = setfield(res, 'time', time);
    
end


if strcmp(flag,'hofcs')
    if nargin<4
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

    z32 = 2^32-1;

    NChannels = floor(sync/head.Resolution);
    tmp   = (0:(NChannels-1))*head.Resolution;
    tcspc = zeros(NChannels,2);
    
    head = pt3Read_t(name1);   
    t_Counts = head.NCounts;
    head = pt3Read_t(name2);        
    t_Counts = t_Counts + head.NCounts;
    
    rate    = [];
    time    = [];
    tst1    = 1; 
    tst2    = 1; 
    cnt     = 0;
    cnt1    = 0; 
    cnt2    = 0; 
 
    max_time  = 0;

    overcount1 = 0;
    overcount2 = 0;
    
    while ((tst1>0)&(tst2>0))

        max_time = max_time + timeframe;

        [head, tmpy1, tmpx1, z1, tst1, overcount1] = pt3Read_t(name1, cnt1+1, overcount1, max_time, CntRate0);
        [head, tmpy2, tmpx2, z2, tst2, overcount2] = pt3Read_t(name2, cnt2+1, overcount2, max_time, CntRate0);

        cnt1 = cnt1 + tst1;
        cnt2 = cnt2 + tst2;

        if ((tst1>0)&(tst2>0))

            tmpy = round([tmpy1+tmpx1; tmpy2+tmpx2]/timeunit);  % in timeunits

            tmpx = [tmpx1; tmpx2];                                  % in ns
            ind = [zeros(length(tmpy1),1); ones(length(tmpy2),1)];

            [tmpy, ord] = sort(tmpy);
            tmpx = tmpx(ord);
            ind = ind(ord);
            tmpy = tmpy-tmpy(1);

            cnt = cnt+length(tmpy);
            dt = (tmpy(end))*timeunit*1e-9;
            rate = [rate length(tmpy)/dt];
            time = [time dt];

            tcspc(:,1) = tcspc(:,1) + mHist(tmpx1, tmp);
            tcspc(:,2) = tcspc(:,2) + mHist(tmpx2, tmp);
            
            tmpy = ceil(tmpy); ind1 = double(ind==0); ind2 = double(ind==1);
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
                        auto(:,j,k,1) = auto(:,j,k,1) + pt32fcs(tmpy, ind1.^j, ind2.^k, Ncasc, Nsub, smooth);
                        auto(:,j,k,2) = auto(:,j,k,2) + pt32fcs(tmpy, ind2.^j, ind1.^k, Ncasc, Nsub, smooth);
                    else
                        auto(:,j,k,1) = auto(:,j,k,1) + tttr2fcs(tmpy, ind1.^j, ind2.^k, Ncasc, Nsub, smooth);
                        auto(:,j,k,2) = auto(:,j,k,2) + tttr2fcs(tmpy, ind2.^j, ind1.^k, Ncasc, Nsub, smooth);
                    end
                end
            end

            subplot('position', [0.875 0.2 0.05 0.6]);
            bar(0,(cnt1+cnt2)/t_Counts);
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
    semilogx(autotime,bld/sum(time)/timeunit*1e3);
    res = struct('autotime', autotime);
    res = setfield(res, 'auto', auto/sum(time)/timeunit*1e3);
    res = setfield(res, 'intensity', intensity/sum(time));
    tmp = 0:head.NChannels;
    res = setfield(res, 'tau', tmp(sum(tcspc')>0)*head.Resolution);
    res = setfield(res, 'tcspc', tcspc(sum(tcspc')>0,:));
    res = setfield(res, 'rate', rate);
    res = setfield(res, 'time', time);
    
end


if strcmp(flag,'flcs')
    if nargin<4
        disp('No master patterns given!');
        return
    end

    Nsub  = 10;
    Ncasc = ceil(log2(3/(sync*1e-9*Nsub)+1));    % 3 sec max correlation time
    autotime = 1e-9*sync*cumsum(reshape(repmat(2.^[0 0:Ncasc-1],Nsub,1),(Ncasc+1)*Nsub,1));
    timeunit = sync;
    
    z32 = 2^32-1;

    head = pt3Read_t(name1);   
    t_Counts = head.NCounts;
    head = pt3Read_t(name2);        
    t_Counts = t_Counts + head.NCounts;
    
    rate    = [];
    time    = [];
    tst1    = 1; 
    tst2    = 1; 
    cnt     = 0;
    cnt1    = 0; 
    cnt2    = 0; 
 
    max_time  = 0;

    overcount1 = 0;
    overcount2 = 0;

    tau = param(:,1);
    param = param(:,2:end);
    
    if nargin<5
        weight = pt3Pro(name1, name2, 'tcspc');
        weight = weight.tcspc(weight.tau>=tau(1) & weight.tau<=tau(end),:);
    end
    semilogy(tau,weight); drawnow;
    
    % filter generation
    % para = [dye1 channel#1; dye1 channel#2; dye2 channel#1; dye2 channel#2];
    for j=1:size(param,2)
        param(:,j) = param(:,j)/sum(param(:,j));
    end % normalization
    master = (param(:,1:2:end)'*diag(1./weight(:,1))*param(:,1:2:end))\(diag(1./weight(:,1))*param(:,1:2:end))';
    master = [master; (param(:,2:2:end)'*diag(1./weight(:,2))*param(:,2:2:end))\(diag(1./weight(:,2))*param(:,2:2:end))']';
    
    auto = zeros(length(autotime),size(master,2)/2,size(master,2)/2);

    while ((tst1>0)&(tst2>0))

        max_time = max_time + timeframe;

        [head, tmpy1, tmpx1, z1, tst1, overcount1] = pt3Read_t(name1, cnt1+1, overcount1, max_time, CntRate0);
        [head, tmpy2, tmpx2, z2, tst2, overcount2] = pt3Read_t(name2, cnt2+1, overcount2, max_time, CntRate0);

        cnt1 = cnt1 + tst1;
        cnt2 = cnt2 + tst2;
        
        if ((tst1>0)&(tst2>0))

            tmpy = round([tmpy1+tmpx1; tmpy2+tmpx2]/timeunit);  % in timeunits

            tmpx = [tmpx1; tmpx2];                                  % in ns

            ind = [zeros(length(tmpy1),1); ones(length(tmpy2),1)];
            
            tst = tmpx>=tau(1) & tmpx<=tau(end);
            tmpy = tmpy(tst);
            tmpx = tmpx(tst);
            ind = ind(tst);

            [tmpy, ord] = sort(tmpy);
            tmpx = tmpx(ord);
            ind = ind(ord);
            tmpy = tmpy-tmpy(1);
            tmpx = round((tmpx-tau(1))/head.Resolution)+1;
            
            cnt = cnt+length(tmpy);
            dt = (tmpy(end))*timeunit*1e-9;
            rate = [rate length(tmpy)/dt];
            time = [time dt];

            for j=1:size(master,2)/2
                for k=1:size(master,2)/2
                    if tmpy(end)>z32
                        auto(:,j,k) = auto(:,j,k) + pt32fcs(tmpy, master(tmpx,j).*double(ind==0), master(tmpx,end/2+k).*double(ind==1), Ncasc, Nsub, smooth);
                        auto(:,j,k) = auto(:,j,k) + pt32fcs(tmpy, master(tmpx,end/2+j).*double(ind==1), master(tmpx,k).*double(ind==0), Ncasc, Nsub, smooth);
                    else
                        auto(:,j,k) = auto(:,j,k) + tttr2fcs(tmpy, master(tmpx,j).*double(ind==0), master(tmpx,end/2+k).*double(ind==1), Ncasc, Nsub, smooth);
                        auto(:,j,k) = auto(:,j,k) + tttr2fcs(tmpy, master(tmpx,end/2+j).*double(ind==1), master(tmpx,k).*double(ind==0), Ncasc, Nsub, smooth);
                    end

                end
            end

            subplot('position', [0.925 0.2 0.025 0.6]);
            bar(0,(cnt1+cnt2)/t_Counts);
            axis([-0.4 0.4 0 1]);
            set(gca,'xtick',[],'ytick',[]);
            subplot('position', [0.1 0.1 0.7 0.8]);
            semilogx(autotime, auto(:,:,1)/sum(time)/timeunit*1e3); 
            drawnow;
        end
    end
    clf 
    semilogx(autotime,auto(:,:,1)/sum(time)/timeunit*1e3);
    res = struct('autotime', autotime);
    res = setfield(res, 'auto', auto/sum(time)/timeunit*1e3);
    res = setfield(res, 'tau', tau);
    res = setfield(res, 'weight', weight);
    res = setfield(res, 'rate', rate);
    res = setfield(res, 'master', master);
    res = setfield(res, 'time', time);
end


if strcmp(flag,'mlcs')
    if nargin<4
        disp('No master patterns given!');
        return
    end

    Nsub  = 10;
    Ncasc = ceil(log2(3/(sync*1e-9*Nsub)+1));    % 3 sec max correlation time
    autotime = 1e-9*sync*cumsum(reshape(repmat(2.^[0 0:Ncasc-1],Nsub,1),(Ncasc+1)*Nsub,1));
    timeunit = sync;
    
    tcspc = zeros(size(param,1),2);
    
    z32 = 2^32-1;

    head = pt3Read_t(name1);   
    t_Counts = head.NCounts;
    head = pt3Read_t(name2);        
    t_Counts = t_Counts + head.NCounts;
    
    rate    = [];
    time    = [];
    tst1    = 1; 
    tst2    = 1; 
    cnt     = 0;
    cnt1    = 0; 
    cnt2    = 0; 
 
    max_time  = 0;

    overcount1 = 0;
    overcount2 = 0;

    tau = param(:,1);
    param = param(:,2:end);
    
    % ML filter generation
    % para = [dye1 channel#1; dye1 channel#2; dye2 channel#1; dye2 channel#2];
    for j=1:size(param,2)
        param(:,j) = param(:,j)/sqrt(sum(param(:,j).^2));
    end % normalization
    master = [param(:,1:2:end) param(:,2:2:end)];
    len = size(master,2)/2;
    
    auto = zeros(length(autotime),len,len);

    while ((tst1>0)&(tst2>0))

        max_time = max_time + timeframe;

        [head, tmpy1, tmpx1, z1, tst1, overcount1] = pt3Read_t(name1, cnt1+1, overcount1, max_time, CntRate0);
        [head, tmpy2, tmpx2, z2, tst2, overcount2] = pt3Read_t(name2, cnt2+1, overcount2, max_time, CntRate0);

        cnt1 = cnt1 + tst1;
        cnt2 = cnt2 + tst2;
        
        if ((tst1>0)&(tst2>0))

            tmpy = round([tmpy1+tmpx1; tmpy2+tmpx2]/timeunit);  % in timeunits

            tmpx = [tmpx1; tmpx2];                                  % in ns

            ind = [zeros(length(tmpy1),1); ones(length(tmpy2),1)];
            
            tst = tmpx>=tau(1) & tmpx<=tau(end);
            tmpy = tmpy(tst);
            tmpx = tmpx(tst);
            ind = ind(tst);

            tcspc(:,1) = tcspc(:,1) + mHist(tmpx(ind==0), tau);
            tcspc(:,2) = tcspc(:,2) + mHist(tmpx(ind==1), tau);

            [tmpy, ord] = sort(tmpy);
            tmpx = tmpx(ord);
            ind = ind(ord);
            tmpy = tmpy-tmpy(1);
            tmpx = round((tmpx-tau(1))/head.Resolution)+1;
            
            cnt = cnt+length(tmpy);
            dt = (tmpy(end))*timeunit*1e-9;
            rate = [rate length(tmpy)/dt];
            time = [time dt];

            for j=1:len
                for k=1:len
                    if tmpy(end)>z32
                        auto(:,j,k) = auto(:,j,k) + pt32fcs(tmpy, master(tmpx,j).*double(ind==0), master(tmpx,end/2+k).*double(ind==1), Ncasc, Nsub, smooth);
                        auto(:,j,k) = auto(:,j,k) + pt32fcs(tmpy, master(tmpx,end/2+k).*double(ind==1), master(tmpx,j).*double(ind==0), Ncasc, Nsub, smooth);
                    else
                        auto(:,j,k) = auto(:,j,k) + tttr2fcs(tmpy, master(tmpx,j).*double(ind==0), master(tmpx,end/2+k).*double(ind==1), Ncasc, Nsub, smooth);
                        auto(:,j,k) = auto(:,j,k) + tttr2fcs(tmpy, master(tmpx,end/2+k).*double(ind==1), master(tmpx,j).*double(ind==0), Ncasc, Nsub, smooth);
                    end
                end
            end

            subplot('position', [0.925 0.2 0.025 0.6]);
            bar(0,(cnt1+cnt2)/t_Counts);
            axis([-0.4 0.4 0 1]);
            set(gca,'xtick',[],'ytick',[]);
            subplot('position', [0.1 0.1 0.7 0.8]);
            semilogx(autotime, auto(:,:,1)/sum(time)/timeunit*1e3); 
            drawnow;
        end
    end
    clf 
    semilogx(autotime,auto(:,:,1)/sum(time)/timeunit*1e3);
    res = struct('autotime', autotime);
    auto = auto/sum(time)/timeunit*1e3;
    res = setfield(res, 'auto0', auto);

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
    res = setfield(res, 'mm', mm);
    auto = permute(auto,[2 3 1]);
    for s=1:size(auto,3)
        auto(:,:,s) = reshape(mm\reshape(squeeze(auto(:,:,s)),len^2,1),len,len);
    end
    auto = permute(auto,[3 1 2]);
    res = setfield(res, 'auto', auto);
    end

    res = setfield(res, 'tau', tau);
    res = setfield(res, 'tcspc', tcspc);
    res = setfield(res, 'rate', rate);
    res = setfield(res, 'master', master);
    res = setfield(res, 'time', time);
end


if strcmp(flag,'antibunch')
    if nargin<4 || isempty(param)
        timeunit = head.Resolution*6.25;
        autotime = (0:2250+1)';
        auto  = zeros(length(autotime),2);
    else
        timeunit = head.Resolution*param(1);
        autotime = (0:param(2)+1)';
        auto  = zeros(length(autotime),2);
    end
    
    if nargin>4
        weight = round(Count0*weight)/Count0*1e9;        
        auto = [auto auto];
    end

    NChannels = floor(sync/head.Resolution);
    tmp   = (0:(NChannels-1))*head.Resolution;
    tcspc = zeros(NChannels,2);
    
    head = pt3Read_t(name1);   
    t_Counts = head.NCounts;
    head = pt3Read_t(name2);        
    t_Counts = t_Counts + head.NCounts;
    
    rate    = [];
    time    = [];
    tst1    = 1; 
    tst2    = 1; 
    cnt     = 0;
    cnt1    = 0; 
    cnt2    = 0; 
 
    max_time  = 0;

    overcount1 = 0;
    overcount2 = 0;
    
    while ((tst1>0)&(tst2>0))

        max_time = max_time + timeframe;

        [head, tmpy1, tmpx1, z1, tst1, overcount1] = pt3Read_t(name1, cnt1+1, overcount1, max_time, CntRate0);
        [head, tmpy2, tmpx2, z2, tst2, overcount2] = pt3Read_t(name2, cnt2+1, overcount2, max_time, CntRate0);

        cnt1 = cnt1 + tst1;
        cnt2 = cnt2 + tst2;

        if ((tst1>0)&(tst2>0))

            tmpy = round([tmpy1+tmpx1; tmpy2+tmpx2]/timeunit);  % in timeunits

            tmpx = [tmpx1; tmpx2];                                  % in ns
            ind = [zeros(length(tmpy1),1); ones(length(tmpy2),1)];

            [tmpy, ord] = sort(tmpy);
            tmpx = tmpx(ord);
            ind = ind(ord);

            cnt = cnt+length(tmpy);
            dt = (tmpy(end) - tmpy(1))*timeunit*1e-9;
            rate = [rate length(tmpy)/dt];
            time = [time dt];

            tcspc(:,1) = tcspc(:,1) + mHist(tmpx1, tmp);
            tcspc(:,2) = tcspc(:,2) + mHist(tmpx2, tmp);

            auto(:,1) = auto(:,1) + mHist(diff(tmpy),autotime,double(ind(1:end-1)==0).*double(ind(2:end)==1));
            auto(:,2) = auto(:,2) + mHist(diff(tmpy),autotime,double(ind(1:end-1)==1).*double(ind(2:end)==0));

            if nargin>4
                tmpy = round([tmpy1+tmpx1; tmpy2+tmpx2-weight]/timeunit);  % in timeunits
                ind = [zeros(length(tmpy1),1); ones(length(tmpy2),1)];
                [tmpy, ord] = sort(tmpy);
                ind = ind(ord);
                auto(:,3) = auto(:,3) + mHist(diff(tmpy),autotime,double(ind(1:end-1)==0).*double(ind(2:end)==1));                
                auto(:,4) = auto(:,4) + mHist(diff(tmpy),autotime,double(ind(1:end-1)==1).*double(ind(2:end)==0));
            end
            
            subplot('position', [0.875 0.2 0.05 0.6]);
            bar(0,(cnt1+cnt2)/t_Counts);
            axis([-0.4 0.4 0 1]);
            set(gca,'xtick',[],'ytick',[]);
            subplot('position', [0.1 0.1 0.7 0.8]);
            plot(autotime(1:end-1)*timeunit,auto(1:end-1,:)/sum(time)); 
            drawnow;
        end

    end
    clf
    
    autotime = autotime(1:end-1)*timeunit;
    auto = auto(1:end-1,:)/sum(time)/timeunit*1e3;
    
    plot(-autotime(end:-1:2),auto(end:-1:2,2),autotime,auto(:,1))
    
    res = struct('autotime', autotime);
    res = setfield(res, 'auto', auto);
    tmp = 0:head.NChannels;
    res = setfield(res, 'tau', tmp(sum(tcspc')>0)*head.Resolution);
    res = setfield(res, 'tcspc', tcspc(sum(tcspc')>0,:));
    res = setfield(res, 'rate', rate);
    res = setfield(res, 'time', time);
    
end


if strcmp(flag,'antibunchm')
    if nargin<4 || isempty(param)
        timeunit = head.Resolution*6.25;
        autotime = (0:2250+1)';
        auto  = zeros(length(autotime),2);
    else
        timeunit = head.Resolution*param(1);
        autotime = (0:param(2)+1)';
        auto  = zeros(length(autotime),2);
    end

    
    if nargin>4
        weight = round(CntRate0*weight)/CntRate0*1e9;
        auto = [auto auto];
    end

    NChannels = floor(sync/head.Resolution);
    tmp   = (0:(NChannels-1))*head.Resolution;
    tcspc = zeros(NChannels,2);
    
    head = pt3Read_t(name1);   
    t_Counts = head.NCounts;
    head = pt3Read_t(name2);        
    t_Counts = t_Counts + head.NCounts;
    
    rate    = [];
    time    = [];
    tst1    = 1; 
    tst2    = 1; 
    cnt     = 0;
    cnt1    = 0; 
    cnt2    = 0; 
 
    max_time  = 0;

    overcount1 = 0;
    overcount2 = 0;

    scnt = 0;
    
    while ((tst1>0)&(tst2>0))

        max_time = max_time + timeframe;

        [head, tmpy1, tmpx1, z1, tst1, overcount1] = pt3Read_t(name1, cnt1+1, overcount1, max_time, CntRate0);
        [head, tmpy2, tmpx2, z2, tst2, overcount2] = pt3Read_t(name2, cnt2+1, overcount2, max_time, CntRate0);

        cnt1 = cnt1 + tst1;
        cnt2 = cnt2 + tst2;

        if ((tst1>0)&(tst2>0))

            scnt = scnt + 1;

            tmpy = round([tmpy1+tmpx1; tmpy2+tmpx2]/timeunit);  % in timeunits

            tmpx = [tmpx1; tmpx2];                                  % in ns
            ind = [zeros(length(tmpy1),1); ones(length(tmpy2),1)];

            [tmpy, ord] = sort(tmpy);
            tmpx = tmpx(ord);
            ind = ind(ord);

            cnt = cnt+length(tmpy);
            dt = (tmpy(end) - tmpy(1))*timeunit*1e-9;
            rate = [rate length(tmpy)/dt];
            time = [time dt];

            tcspc(:,1,scnt) = mHist(tmpx1, tmp);
            tcspc(:,2,scnt) = mHist(tmpx2, tmp);

            auto(:,1,scnt) = mHist(diff(tmpy),autotime,double(ind(1:end-1)==0).*double(ind(2:end)==1));
            auto(:,2,scnt) = mHist(diff(tmpy),autotime,double(ind(1:end-1)==1).*double(ind(2:end)==0));

            if nargin>4
                tmpy = round([tmpy1+tmpx1; tmpy2+tmpx2-weight]/timeunit);  % in timeunits
                ind = [zeros(length(tmpy1),1); ones(length(tmpy2),1)];
                [tmpy, ord] = sort(tmpy);
                ind = ind(ord);
                auto(:,3,scnt) = mHist(diff(tmpy),autotime,double(ind(1:end-1)==0).*double(ind(2:end)==1));                
                auto(:,4,scnt) = mHist(diff(tmpy),autotime,double(ind(1:end-1)==1).*double(ind(2:end)==0));
            end
            
            subplot('position', [0.875 0.2 0.05 0.6]);
            bar(0,(cnt1+cnt2)/t_Counts);
            axis([-0.4 0.4 0 1]);
            set(gca,'xtick',[],'ytick',[]);
            subplot('position', [0.1 0.1 0.7 0.8]);
            plot(autotime(1:end-1)*timeunit,sum(auto(1:end-1,:,:),3)/sum(time)); 
            drawnow;
        end

    end
    clf
    
    autotime = autotime(1:end-1)*timeunit;
    for j=1:size(auto,3) 
        auto(:,:,j) = auto(:,:,j)/time(j)/timeunit*1e3;
    end
    auto = auto(1:end-1,:,:);
    
    plot(-autotime(end:-1:2),mean(auto(end:-1:2,2,:),3),autotime,mean(auto(:,1,:),3))
    
    res = struct('autotime', autotime);
    res = setfield(res, 'auto', auto);
    tmp = 0:head.NChannels;
    res = setfield(res, 'tau', tmp*head.Resolution);
    res = setfield(res, 'tcspc', tcspc);
    res = setfield(res, 'rate', rate);
    res = setfield(res, 'time', time);
    
end


if strcmp(flag,'fimda')
    Ncasc = 10;
    Nsub = 10;        
    nums = 10000;
	timeunit = sync;

    autotime = 1e-9*timeunit*cumsum(reshape(repmat(2.^(0:Ncasc-1),Nsub,1),Ncasc*Nsub,1));
    fimda = zeros(nums+1,length(autotime));
    tcspc = zeros(1,head.NChannels);

    NChannels = floor(sync/head.Resolution);
    tmp   = (0:(NChannels-1))*head.Resolution;
    tcspc = zeros(NChannels,2);
    
    head = pt3Read_t(name1);   
    t_Counts = head.NCounts;
    head = pt3Read_t(name2);        
    t_Counts = t_Counts + head.NCounts;
    
    rate    = [];
    time    = [];
    tst1    = 1; 
    tst2    = 1; 
    cnt     = 0;
    cnt1    = 0; 
    cnt2    = 0; 
 
    max_time  = 0;

    overcount1 = 0;
    overcount2 = 0;
    
    while ((tst1>0)&(tst2>0))

        max_time = max_time + timeframe;

        [head, tmpy1, tmpx1, z1, tst1, overcount1] = pt3Read_t(name1, cnt1+1, overcount1, max_time, CntRate0);
        [head, tmpy2, tmpx2, z2, tst2, overcount2] = pt3Read_t(name2, cnt2+1, overcount2, max_time, CntRate0);

        cnt1 = cnt1 + tst1;
        cnt2 = cnt2 + tst2;

        if ((tst1>0)&(tst2>0))

            tmpy = round([tmpy1+tmpx1; tmpy2+tmpx2]/timeunit);  % in timeunits

            tmpx = [tmpx1; tmpx2];                                  % in ns
            ind = [zeros(length(tmpy1),1); ones(length(tmpy2),1)];

            [tmpy, ord] = sort(tmpy);
            tmpx = tmpx(ord);
            ind = ind(ord);

            cnt = cnt+length(tmpy);
            dt = (tmpy(end) - tmpy(1))*timeunit;
            rate = [rate length(tmpy)/dt];
            time = [time dt];

            tcspc(:,1) = tcspc(:,1) + mHist(tmpx1, tmp);
            tcspc(:,2) = tcspc(:,2) + mHist(tmpx2, tmp);

            fimda = fimda + tttr2fimda(tmpy,ones(size(tmpy)),Ncasc,Nsub,nums);
            subplot('position', [0.875 0.2 0.05 0.6]);
            bar(0,(cnt1+cnt2)/t_Counts);
            axis([-0.4 0.4 0 1]);
            set(gca,'xtick',[],'ytick',[]);
            subplot('position', [0.1 0.1 0.7 0.8]);
            loglog(0:nums,fimda); 
            drawnow;
        end

    end
    clf
    
    loglog(0:nums,fimda); 
    
    res = struct('nums', 0:nums);
    res = setfield(res, 'fimda', fimda);
    res = setfield(res, 'autotime', autotime);
    tmp = 0:head.NChannels;
    res = setfield(res, 'tau', tmp(sum(tcspc')>0)*head.Resolution);
    res = setfield(res, 'tcspc', tcspc(sum(tcspc')>0,:));
    res = setfield(res, 'rate', rate);
    res = setfield(res, 'time', time);
    
end


head.NCounts = cnt;
head.CntRate0 = CntRate0;    
warning on MATLAB:DeprecatedLogicalAPI