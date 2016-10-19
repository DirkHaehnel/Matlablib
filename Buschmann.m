function [res, head] = Buschmann(name, flag, param, weight);

warning off MATLAB:DeprecatedLogicalAPI

timeframe = 6e10; % time-frame to read photons from file in ns
smooth = 1; % smoothed autocorrelation
close all   

head  = pt3ReadBuschmann(name);
CntRate0 = head.CntRate0;

sync = 1e9/CntRate0; 


if strcmp(flag,'tcspc')
    NChannels = floor(sync/head.Resolution);
    tau   = (0:(NChannels-1))*head.Resolution;
    tcspc = zeros(NChannels,1);

	timeunit = sync;
    t_Counts = head.NCounts;

    rate    = [];
    time    = [];
    tst    = 1;
    cnt     = 0;

    max_time  = 0;

    overcount = 0;

    while (tst>0)

        max_time = max_time + timeframe;

        [head, tmpy, tmpx, ind, tst, overcount] = pt3ReadBuschmann(name, cnt+1, overcount, max_time, CntRate0);

        cnt = cnt + tst;

        if (tst>0)

            tmpy = round((tmpy+tmpx)/timeunit);  % in timeunits
            tmpy = tmpy-tmpy(1);
            dt = (tmpy(end))*timeunit*1e-9;
            rate = [rate length(tmpy)/dt];
            time = [time dt];

            tcspc = tcspc + mHist(tmpx(ind==2), tau);

            subplot('position', [0.875 0.2 0.05 0.6]);
            bar(0,cnt/t_Counts);
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
    res = setfield(res, 'rate', rate);
    res = setfield(res, 'time', time);
    
end


if strcmp(flag,'xfcs')
    if nargin<4
        Nsub  = 20;
        Ncasc = ceil(log2(1/(sync*1e-9*Nsub)+1));    % 3 sec max correlation time
        autotime = 1e-9*sync*cumsum(reshape(repmat(2.^[0 0:Ncasc-1],Nsub,1),(Ncasc+1)*Nsub,1));
        auto  = zeros(length(autotime),1);
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
        auto  = zeros(length(autotime),1);
    end
    
    z32 = 2^32-1;

    NChannels = floor(sync/head.Resolution);
    tau   = (0:(NChannels-1))*head.Resolution;
    tcspc = zeros(NChannels,2);
    
    head = pt3Read_t(name);   
    t_Counts = head.NCounts;
    
    rate    = [];
    time    = [];
    tst    = 1; 
    cnt     = 0;
 
    max_time  = 0;

    overcount = 0;
    
    while (tst>0)

        max_time = max_time + timeframe;

        [head, tmpy, tmpx, ind, tst, overcount] = pt3ReadBuschmann(name, cnt+1, overcount, max_time, CntRate0);

        cnt = cnt + tst;

        if (tst>0)

            tmpy = round((tmpy+tmpx)/timeunit);  % in timeunits
            tmpy = tmpy-tmpy(1);

            dt = (tmpy(end))*timeunit*1e-9;
            rate = [rate length(tmpy)/dt];
            time = [time dt];

           tcspc(:,1) = tcspc(:,1) + mHist(tmpx(ind==2), tau);

           if tmpy(end)>z32
               ind32 = tmpy<=z32;
               auto(:,1) = auto(:,1) + tttr2fcs(tmpy(ind32), double((ind(ind32)==2).^2), double((ind(ind32)==2).^2), Ncasc, Nsub, smooth);
           else
                auto(:,1) = auto(:,1) + tttr2fcs(tmpy, double(ind==2), double(ind==2), Ncasc, Nsub, smooth);
           end

            subplot('position', [0.875 0.2 0.05 0.6]);
            bar(0,cnt/t_Counts);
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


if strcmp(flag,'flcs')
    if nargin<3
        disp('No master patterns given!');
        return
    end

    Nsub  = 10;
    Ncasc = ceil(log2(1/(sync*1e-9*Nsub)+1));    % 1 sec max correlation time
    autotime = 1e-9*sync*cumsum(reshape(repmat(2.^[0 0:Ncasc-1],Nsub,1),(Ncasc+1)*Nsub,1));
    timeunit = sync;
    
    z32 = 2^32-1;

    head = pt3Read_t(name);   
    t_Counts = head.NCounts;
    
    rate    = [];
    time    = [];
    tst    = 1; 
    cnt     = 0;
 
    max_time  = 0;

    overcount = 0;

    tau = param(:,1);
    param = param(:,2:end);
    
    if nargin<4
        weight = Buschmann(name, 'tcspc');
        weight = weight.tcspc(weight.tau>=tau(1) & weight.tau<=tau(end),:);
    end
    semilogy(tau,weight); drawnow;
    
    % filter generation
    % para = dye1 dye2;
    for j=1:size(param,2)
        param(:,j) = param(:,j)/sum(param(:,j));
    end % normalization
    master = ((param'*diag(1./weight)*param)\(diag(1./weight)*param)')';
    
    auto = zeros(length(autotime),size(master,2),size(master,2));

    while (tst>0)

        max_time = max_time + timeframe;

        [head, tmpy, tmpx, ind, tst, overcount] = pt3ReadBuschmann(name, cnt+1, overcount, max_time, CntRate0);

        cnt = cnt + tst;
       
        if (tst>0)

            tmpy = round((tmpy+tmpx)/timeunit);  % in timeunits

            tst = tmpx>=tau(1) & tmpx<=tau(end);
            tmpy = tmpy(tst);
            tmpx = tmpx(tst);
            ind = ind(tst);

            dt = (tmpy(end))*timeunit*1e-9;
            rate = [rate length(tmpy)/dt];
            time = [time dt];

            tmpx = round((tmpx-tau(1))/head.Resolution)+1;
            for j=1:size(master,2)
                for k=1:size(master,2)
                    if tmpy(end)>z32
                        ind32 = tmpy<=z32;
                        auto(:,j,k) = auto(:,j,k) + tttr2fcs(tmpy(ind32), master(tmpx(ind32),j), master(tmpx(ind32),k), Ncasc, Nsub, smooth);
                    else
                        auto(:,j,k) = auto(:,j,k) + tttr2fcs(tmpy, master(tmpx,j), master(tmpx,k), Ncasc, Nsub, smooth);
                    end
                end
            end

            subplot('position', [0.925 0.2 0.025 0.6]);
            bar(0,cnt/t_Counts);
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
