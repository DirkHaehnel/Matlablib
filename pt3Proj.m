function [res, head] = pt3Proj(name1, name2, flag, para);

warning off MATLAB:DeprecatedLogicalAPI

timewin = 1; % time window of data processing in s
timeframe = 100e9; % time-frame to read photons from file in ns
smooth = 1; % smoothed autocorrelation
close all

head = pt3Read_t(name1);
timeunit = head.Resolution; % in ns

if strcmp(flag,'tcspc')
    tmp   = 0:head.Resolution:(head.NChannels-1)*head.Resolution;
    tcspc = zeros(head.NChannels,2);

    max_time  = 0;
    tst1    = 1;
    tst2    = 1;
    cnt1    = 0;
    cnt2    = 0;
    overcount1 = 0;
    overcount2 = 0;

    while ((tst1>0)&(tst2>0))

        max_time = max_time + timeframe;

        [head, tmpy1, tmpx1, tst1, overcount1] = pt3Read_t(name1, cnt1+1, overcount1, max_time);
        [head, tmpy2, tmpx2, tst2, overcount2] = pt3Read_t(name2, cnt2+1, overcount2, max_time);

        cnt1 = cnt1 + tst1;
        cnt2 = cnt2 + tst2;

        if ((tst1>0)&(tst2>0))

            tcspc(:,1) = tcspc(:,1) + mHist(tmpx1, tmp);
            tcspc(:,2) = tcspc(:,2) + mHist(tmpx2, tmp);

            subplot('position', [0.875 0.2 0.05 0.6]);
            bar(0,cnt2/head.NCounts);
            axis([-0.4 0.4 0 1]);
            set(gca,'xtick',[],'ytick',[]);
            subplot('position', [0.1 0.1 0.7 0.8]);
            semilogy(tmp,tcspc);
            drawnow;
        end
    end
    clf
    cnt = cnt1+cnt2;
    res = struct('tau', tmp);
    res = setfield(res, 'tcspc', tcspc);
    semilogy(res.tau,res.tcspc);
end


if strcmp(flag,'xfcs')
    %Nsub  = 10;
    %Ncasc = ceil(log2(1/(timeunit*1e-9*Nsub)+1));    % 1 sec max correlation time
    %autotime = 1e-9*timeunit*cumsum(reshape(repmat(2.^[0 0:Ncasc-1],Nsub,1),(Ncasc+1)*Nsub,1));
    Nsub = 100;
    Ncasc = 10;
    autotime = 1e-9*timeunit*cumsum(reshape(repmat(2.^(0:Ncasc-1),Nsub,1),Ncasc*Nsub,1));
    auto  = zeros(length(autotime),2);

    tmp   = 0:head.Resolution:(head.NChannels-1)*head.Resolution;
    tcspc = zeros(head.NChannels,2);

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

        [head, tmpy1, tmpx1, tst1, overcount1] = pt3Read_t(name1, cnt1+1, overcount1, max_time);
        [head, tmpy2, tmpx2, tst2, overcount2] = pt3Read_t(name2, cnt2+1, overcount2, max_time);

        cnt1 = cnt1 + tst1;
        cnt2 = cnt2 + tst2;

        if ((tst1>0)&(tst2>0))

            tmpy = round([tmpy1+tmpx1; tmpy2+tmpx2]/timeunit);          % in timeunits

            tmpx = [tmpx1; tmpx2];                                      % in ns
            ind  = [zeros(length(tmpy1),1); ones(length(tmpy2),1)];

            [tmpy, ord]  = sort(tmpy);
            ind          = ind(ord);

            cnt  = cnt+length(tmpy);
            dt   = tmpy(end) - tmpy(1);
            rate = [rate length(tmpy)/dt];
            time = [time dt];

            tcspc(:,1) = tcspc(:,1) + mHist(tmpx1, tmp);
            tcspc(:,2) = tcspc(:,2) + mHist(tmpx2, tmp);

%             auto(:,1) = auto(:,1) + pt32fcs(tmpy, double(ind==0), double(ind==1), Ncasc, Nsub, smooth);
%             auto(:,2) = auto(:,2) + pt32fcs(tmpy, double(ind==1), double(ind==0), Ncasc, Nsub, smooth);
            auto = tttr2xfcs(tmpy, double(ind==0), double(ind==1), Ncasc, Nsub);
            auto = tttr2xfcs(tmpy, double(ind==1), double(ind==0), Ncasc, Nsub);

            subplot('position', [0.875 0.2 0.05 0.6]);
            bar(0,(cnt1+cnt2)/t_Counts);
            axis([-0.4 0.4 0 1]);
            set(gca,'xtick',[],'ytick',[]);
            subplot('position', [0.1 0.1 0.7 0.8]);
            semilogx(autotime,auto/sum(time)/timeunit/1e-9); %, autotime, 0.5*max(max(auto/sum(time)))*(diff([1; diff([0; autotime])])>0))
            drawnow;
        end
        break
    end
    clf

    semilogx(autotime,auto/sum(time)/timeunit/1e-9);
    res = struct('autotime', autotime);
    res = setfield(res, 'auto', auto/sum(time)/timeunit/1e-9);
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

    tau1 = para(:,1);
    tau2 = para(:,2);
    para = para(:,3:end);

    tcspc = pt3Proj(name1,name2,'tcspc');
    bin = tcspc.tau';
    tcspc = tcspc.tcspc;
    semilogy(bin,tcspc); drawnow;

    % filter generation
    % para = [dye1 channel#1; dye1 channel#2; dye2 channel#1; dye2 channel#2];
    for j=1:size(para,2)
        para(:,j) = para(:,j)/sum(para(:,j));
    end % normalization
    weight = [tcspc(tau1(1)<=bin & bin<=tau1(end),1) tcspc(tau2(1)<=bin & bin<=tau2(end),2)];
    master = (para(:,1:2:end)'*diag(1./weight(:,1))*para(:,1:2:end))\(diag(1./weight(:,1))*para(:,1:2:end))';
    master = [master; (para(:,2:2:end)'*diag(1./weight(:,2))*para(:,2:2:end))\(diag(1./weight(:,2))*para(:,2:2:end))']';
    dtau = mean(diff(tau1));

    Nsub  = 10;
    Ncasc = ceil(log2(1/(timeunit*1e-9*Nsub)+1));    % 1 sec max correlation time
    autotime = 1e-9*timeunit*cumsum(reshape(repmat(2.^[0 0:Ncasc-1],Nsub,1),(Ncasc+1)*Nsub,1));
    auto  = zeros(length(autotime),size(master,2)/2,size(master,2)/2);

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

        [head, tmpy1, tmpx1, tst1, overcount1] = pt3Read_t(name1, cnt1+1, overcount1, max_time);
        [head, tmpy2, tmpx2, tst2, overcount2] = pt3Read_t(name2, cnt2+1, overcount2, max_time);

        cnt1 = cnt1 + tst1;
        cnt2 = cnt2 + tst2;

        if ((tst1>0)&(tst2>0))

            tmpy = round([tmpy1+tmpx1; tmpy2+tmpx2]/timeunit);          % in timeunits

            tmpx = [tmpx1; tmpx2];                                      % in ns
            ind  = [zeros(length(tmpy1),1); ones(length(tmpy2),1)];

            [tmpy, ord]  = sort(tmpy);
            ind          = ind(ord);

            cnt  = cnt+length(tmpy);
            dt   = tmpy(end) - tmpy(1);
            rate = [rate length(tmpy)/dt];
            time = [time dt];

            tmp = (ind==0 & (tmpx<tau1(1) | tau1(end)<tmpx)) | (ind==1 & (tmpx<tau2(1) | tau2(end)<tmpx));
            tmpy(tmp) = [];
            tmpx(tmp) = [];
            ind(tmp) = [];
            tmpx(ind==0) = tmpx(ind==0)-tau1(1) + dtau;
            tmpx(ind==1) = tmpx(ind==1)-tau2(1) + dtau;
            tmpx = round(tmpx/dtau);
            ind = double(ind);
            for j=1:size(master,2)/2
                for k=1:size(master,2)/2
                    auto(:,j,k) = auto(:,j,k) + pt32fcs(tmpy, master(tmpx,j).*ind, ...
                        master(tmpx,end/2+k).*(1-ind), Ncasc, Nsub, smooth);
                end
            end

            subplot('position', [0.875 0.2 0.05 0.6]);
            bar(0,(cnt1+cnt2)/t_Counts);
            axis([-0.4 0.4 0 1]);
            set(gca,'xtick',[],'ytick',[]);
            subplot('position', [0.1 0.1 0.7 0.8]);
            semilogx(autotime,auto(:,:,1)/sum(time)/timeunit/1e-9); 
            drawnow;
        end
    end
    clf

    semilogx(autotime,auto(:,:,1)/sum(time)/timeunit/1e-9); 
    res = struct('autotime', autotime);
    res = setfield(res, 'auto', auto/sum(time)/timeunit/1e-9);
    res = setfield(res, 'tau', bin');
    res = setfield(res, 'tcspc', tcspc);
    res = setfield(res, 'master', master);
    res = setfield(res, 'weight', weight);
end


if strcmp(flag,'antibunch')

    Nsub  = 1000;
    Ncasc = 1;           % 100ns max correlation time

    autotime = 1E-10*cumsum(reshape(repmat(2.^[0 0:Ncasc-1],Nsub,1),(Ncasc+1)*Nsub,1))-1E-7;
    auto     = zeros(length(autotime),2);

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

    timeframe = 100E9;   % time-frame to read photons from file in ns
    max_time  = 0;

    overcount1 = 0;
    overcount2 = 0;

    while ((tst1>0)&(tst2>0))

        max_time = max_time + timeframe;

        [head, tmpy1, tmpx1, tst1, overcount1] = pt3Read_t(name1, cnt1+1, overcount1, max_time);
        [head, tmpy2, tmpx2, tst2, overcount2] = pt3Read_t(name2, cnt2+1, overcount2, max_time);

        cnt1 = cnt1 + tst1;
        cnt2 = cnt2 + tst2;

        if ((tst1>0)&(tst2>0))

            tmpy  = round(10*[tmpy1+tmpx1; tmpy2+tmpx2]);              % in 100 ps
            cnt   = cnt+length(tmpy);
            dt    = tmpy(end) - tmpy(1);
            rate  = [rate length(tmpy)/dt];
            time  = [time dt];

            tmps1 = round(10*[100+tmpy1+tmpx1; tmpy2+tmpx2]);          % in 100 ps
            tmps2 = round(10*[tmpy1+tmpx1; 100+tmpy2+tmpx2]);          % in 100 ps
            ind1  = [zeros(length(tmpy1),1); ones(length(tmpy2),1)];
            ind2  = [zeros(length(tmpy1),1); ones(length(tmpy2),1)];

            [tmps1, ord1]  = sort(tmps1);
            [tmps2, ord2]  = sort(tmps2);
            ind1           = ind1(ord1);
            ind2           = ind2(ord2);

            tc = tmps1; auto(:,1) = auto(:,1) + pt32fcs(tc, double(ind1==1), double(ind1==0), Ncasc, Nsub, smooth);
            tc = tmps2; auto(:,2) = auto(:,2) + pt32fcs(tc, double(ind2==0), double(ind2==1), Ncasc, Nsub, smooth);

            subplot('position', [0.875 0.2 0.05 0.6]);
            bar(0,(cnt1+cnt2)/t_Counts);
            axis([-0.4 0.4 0 1]);
            set(gca,'xtick',[],'ytick',[]);
            subplot('position', [0.1 0.1 0.7 0.8]);
            plot(autotime,1E10*auto/sum(time)); %, autotime, 0.5*max(max(auto/sum(time)))*(diff([1; diff([0; autotime])])>0))
            drawnow;
        end
    end
    clf

    auto = 1E10*auto/sum(time);
    auto(:,2) = auto(end:-1:1,2);

    yy1 = smooth(mean(auto,2),20,'rloess');

    plot(autotime,[auto yy1]);

    res = struct('autotime', autotime);
    res = setfield(res, 'auto', auto);
    res = setfield(res, 'smooth', yy1);
    res = setfield(res, 'rate', rate);
    res = setfield(res, 'time', time);
end


head.NCounts = cnt;

warning on MATLAB:DeprecatedLogicalAPI