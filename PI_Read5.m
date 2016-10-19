function [tag, tim, tcspc, bin] = PI_Read(name, flag, filter, numflags)

[tmp, tmp, tmp, head] = tttrRead(name,[1 1]);

para = [head.ScanTimePerPix*2e3, head.ScanWidthX, head.ScanWidthY];

if nargin<4
    numflags = 2;
end

if nargin<3 || isempty(filter)
    filt = [];
else
    filt = sort(filter(1,:)); 
end

maxCounts = 1e5;

nch  = 2;
bin  = 1:4096;

if strcmp(flag, '2focus')
    
    tst   = 1;
    cnt   = 0;
    tcspc = zeros(head.NChannels,2);
    tmp   = 1:head.NChannels;    

    h = waitbar(0,'Please wait ...'); drawnow
    
    while tst>0
        [tmpy, flv, tmpx, head, tst, overcount, marker] = tttrRead(name, [cnt+1 maxCounts]);        
        
        flv(marker>0)  = [];
        tmpx(marker>0) = [];        
        
        if tst>0
            cnt = cnt+tst;
            tcspc(:,1) = tcspc(:,1) + mHist(tmpx(flv==0), tmp);
            tcspc(:,2) = tcspc(:,2) + mHist(tmpx(flv==1), tmp);            
            waitbar(cnt/head.NCounts);
        end
    end
    close(h);

% Assign windows in TCSPC-histgramm to each laser/detector
% Windows are constructed so, that the maximum is in the 30th channel of
% each windwow. The windows correspond to the following assignment:
% 1st window: laser 1, SPAD 1
% 2nd window: laser 1, SPAD 2
% 3rd window: laser 2, SPAD 1
% 4th window: laser 2, SPAD 2

    ind = sum(tcspc,2)==0;
    
    tmp(ind) = [];
    tcspc(ind,:) = [];
    tmp(1:30)=[];
    tcspc(1:30,:) = [];

    tmpa  = tcspc(:,1);
    tmpb  = tcspc(:,2);
    bin   = tmp';

    win   = (1:round(length(bin)/3.5))-30;

    tcspc = zeros(length(win),4); 
    tau   = tcspc;
    
    for j=1:2
        [pos pos] = max(tmpa);
        ind = pos + win;
        if ind(1)<=0
            win(1:sum(ind<=0))     = [];
            tcspc(1:sum(ind<=0),:) = [];
            tau(1:sum(ind<=0),:)   = [];
            ind = pos + win;
        end
        if ind(end)>length(tmpa)
            win(end-sum(ind>length(tmpa)):end)     = [];
            tcspc(end-sum(ind>length(tmpa)):end,:) = [];
            tau(end-sum(ind>length(tmpa)):end,:)   = [];
            ind = pos + win;
        end
        tcspc(:,2*j-1) = tmpa(ind);
        tau(:,2*j-1)   = bin(ind);
        tmpa(ind)      = 0;
        
        [pos pos] = max(tmpb);
        ind = pos + win;
        if ind(1)<=0
            win(1:sum(ind<=0))     = [];
            tcspc(1:sum(ind<=0),:) = [];
            tau(1:sum(ind<=0),:)   = [];
            ind = pos + win;
        end
        if ind(end)>length(tmpb)
            win(end-sum(ind>length(tmpb)):end)     = [];
            tcspc(end-sum(ind>length(tmpb)):end,:) = [];
            tau(end-sum(ind>length(tmpb)):end,:)   = [];
            ind = pos + win;
        end
        tcspc(:,2*j) = tmpb(ind);
        tau(:,2*j)   = bin(ind);
        tmpb(ind)    = 0;
    end
    
    [ind, ind] = sort(tau(1,1:2:end));
    ind        = 2*(ind-1)+1;
    tau(:,1:2:end)   = tau(:,ind);
    tcspc(:,1:2:end) = tcspc(:,ind);
    
    [ind,ind] = sort(tau(1,2:2:end));
    ind       = 2*ind;
    tau(:,2:2:end)   = tau(:,ind);
    tcspc(:,2:2:end) = tcspc(:,ind);
    
    if ~isempty(filt)
        filter = round(filt/head.Resolution) + 30;  % convert to channels
        if filter(1)<1
            filter(1) = 1;
        end;
        if filter(2)<2
            filter(2) = 2;
        end;
        if filter(1)>size(tau,1)-1
            filter(1) = size(tau,1)-1;
        end;
        if filter(2)>size(tau,1)-1
            filter(2) = size(tau,1);
        end;
        ind = filter(1):filter(2);
        tau = tau(ind,:);
        tcspc = tcspc(ind);
    end;

    for j=1:4                                       
        tcspc(:,j) = tcspc(:,j)-min(tcspc(:,j));    
    end                                                                                                 
    
    h = waitbar(0,'Please wait ...'); drawnow
        
    cnts      = 0;
    starttime = 0;
    line      = 0;
    tmpy = [];
    tmpf = [];
    tmpt = [];
    tmpm = [];
    tag  = [];
    tim  = [];
    
    while cnts < head.NCounts
        [y, flag, tcspc, head, num, overcount, marker] = tttrRead(name, [cnts+1 maxCounts]);

        y = y + starttime;
        starttime = starttime + overcount;

        marker = bitget(marker,3);

        cnts = cnts + num;
        waitbar(cnts/head.NCounts);

        tmpy = [tmpy; y];
        tmpf = [tmpf; flag];
        tmpt = [tmpt; tcspc];
        tmpm = [tmpm; marker];

        tmp = find(tmpm~=0);
        while numel(tmp)>=numflags
            line       = line+1
            
            p          = (1:para(3)).*head.ScanResolution;
            line_start = tmpy(tmp(1));
            line_t     = tmpy(tmp(1:numflags)) - line_start;
            line_p     = (0:1/(numflags-1):1)'.*para(3).*head.ScanResolution;
            
            tc = polyfit(line_p, line_t, min([3,numflags-1]));
            t  = polyval(tc, p); 
%           t = interp1(line_p, line_t, p, 'cubic');
            
            rsty      = tmpy(tmp(numflags)+1:end);
            rstf      = tmpf(tmp(numflags)+1:end);
            rsttcspc  = tmpt(tmp(numflags)+1:end);
            rstmarker = tmpm(tmp(numflags)+1:end);
            
            tmpy = [tmpy(tmp(1):tmp(numflags))] - tmpy(tmp(1));
            tmpt = [tmpt(tmp(1):tmp(numflags))];
            tmpf = [tmpf(tmp(1):tmp(numflags))];
            tmpm = [tmpm(tmp(1):tmp(numflags))];
            
            tmpy(tmpm~=0) = [];
            tmpt(tmpm~=0) = [];
            tmpf(tmpm~=0) = [];
            tmpm(tmpm~=0) = [];


            tmptag = zeros(para(3),2);
            tmptim = zeros(para(3),2);

            for n = 1:para(3)
                tmp = numel(tmpy(tmpy<t(n)));
                
                if tmp >0
                    ttch  = tmpt(1:tmp);
                    ttmpf = tmpf(1:tmp);

                    ind = zeros(tmp,4);
                    ind(:,1) = (ttmpf==0 & ttch>=tau(1,1) & tau(end,1)>=ttch);  % photons of laser 1 detected by SPAD 1
                    ind(:,2) = (ttmpf==1 & ttch>=tau(1,2) & tau(end,2)>=ttch);  % photons of laser 1 detected by SPAD 2
                    ind(:,3) = (ttmpf==0 & ttch>=tau(1,3) & tau(end,3)>=ttch);  % photons of laser 2 detected by SPAD 1
                    ind(:,3) = (ttmpf==1 & ttch>=tau(1,4) & tau(end,4)>=ttch);  % photons of laser 2 detected by SPAD 2

                    tmptag(n,1) = sum(ind(:,1))+sum(ind(:,2));
                    tmptag(n,2) = sum(ind(:,3))+sum(ind(:,4));

                    ttch(ind(:,1)==1) = ttch(ind(:,1)==1) - tau(30,1);
                    ttch(ind(:,2)==1) = ttch(ind(:,2)==1) - tau(30,2);
                    ttch(ind(:,3)==1) = ttch(ind(:,3)==1) - tau(30,3);
                    ttch(ind(:,4)==1) = ttch(ind(:,4)==1) - tau(30,4);
                    ttch(ttch<0) = 0;
                    ttch = ttch*head.Resolution; % in ns

                    tmptim(n,1) = (sum(ttch(ind(:,1)==1))+sum(ttch(ind(:,2)==1)));
                    tmptim(n,2) = (sum(ttch(ind(:,3)==1))+sum(ttch(ind(:,4)==1)));

                    tmpy(1:tmp)  = [];
                    tmpt(1:tmp)  = [];
                    tmpf(1:tmp)  = [];
                    tmpm(1:tmp)  = [];
                    ind(1:tmp,:) = [];
                end
            end;
            
            tag = [tag; tmptag];
            tim = [tim; tmptim./tmptag];

            tmpy = [rsty];
            tmpf = [rstf];
            tmpt = [rsttcspc];
            tmpm = [rstmarker];

            rsty      = [];
            rstf      = [];
            rsttcspc  = [];
            rstmarker = [];

            tmp = find(tmpm~=0);
        end;
    end

    close(h);

    tag = reshape(tag, para(2), para(3), 2);
    tag(:, 2:2:end, :)  = tag(end:-1:1, 2:2:end, :);
    tag = permute(tag, [2 1 3]);

    nx        = head.ScanWidthX;
    ny        = head.ScanWidthY;

    if ~isempty(tim)&&(size(tim,1)==nx*ny)
        tim = reshape(tim, para(2), para(3), numel(tim)/prod(para(2:3))/2, 2);
        tim(:, 2:2:end, :, :) = tim(end:-1:1, 2:2:end, :, :);
        tim = permute(tim, [2 1 3:ndims(tim)]);
        tim = squeeze(tim);
    end;

    figure;
    x = head.ScanStartX+(1:size(tag,1))*head.ScanResolution;
    y = head.ScanStartY+(1:size(tag,2))*head.ScanResolution;
    imagesc(y,x,tag(:,:,1));
    set(gca,'DataAspectRatio', [1,1,1]);
    set(gca,'XDir','reverse');
    set(gca,'YDir','reverse');
    xlabel('x / µm');
    ylabel('y / µm');
    figure;
    x = head.ScanStartX+(1:size(tag,1))*head.ScanResolution;
    y = head.ScanStartY+(1:size(tag,2))*head.ScanResolution;
    imagesc(y,x,tag(:,:,2));
    set(gca,'DataAspectRatio', [1,1,1]);
    set(gca,'XDir','reverse');
    set(gca,'YDir','reverse');
    xlabel('x / µm');
    ylabel('y / µm');

end


if strcmp(flag, 'channel')    
    cnts      = 0;
    in_line   = false;
    starttime = 0;

    tmpy = [];
    tmpf = [];
    tmpt = [];
    tmpm = [];
    tag  = [];
    tim  = [];


    while cnts < head.NCounts
        [y, flag, tcspc, head, num, overcount, marker] = tttrRead(name, [cnts+1 maxCounts]);

        y = y + starttime;
        starttime = starttime + overcount;

        marker = bitget(marker,3);

        cnts = cnts + num;
        waitbar(cnts/head.NCounts);

        tmpy = [tmpy; y];
        tmpf = [tmpf; flag];
        tmpt = [tmpt; tcspc];
        tmpm = [tmpm; marker];

        tmp = find(tmpm~=0);
        while numel(tmp)>0
            if ~in_line
                tmpy(1:tmp(1)) = [];
                tmpf(1:tmp(1)) = [];
                tmpt(1:tmp(1)) = [];
                tmpm(1:tmp(1)) = [];
                in_line        = true;
            else
                rsty      = tmpy(tmp(1)+1:end);
                rstf      = tmpf(tmp(1)+1:end);
                rsttcspc  = tmpt(tmp(1)+1:end);
                rstmarker = tmpm(tmp(1)+1:end);

                tmpy = [tmpy(1:tmp(1)-1) - tmpy(1)];
                tmpt = [tmpt(1:tmp(1)-1)];
                tmpf = [tmpf(1:tmp(1)-1)];
                tmpm = [tmpm(1:tmp(1)-1)];

                tmptag = tttr2bin(tmpy, 1+tmpy(end)/para(3), tmpf);
                tag    = [tag; tmptag];

                for jch=1:nch
                    tcs = tmpt(tmpf==jch-1);
                    if nargout>2
                        tch(:,jch) = tch(:,jch) + mhist(tcs,bin);
                    end
                    tmp = cumsum([1; tmptag(:,jch)]);
                    if nargin<2
                        tcs = cumsum([0; tcs]);
                        tmp = cumsum([1; tmptag(:,jch)]);
                        tim = [tim; (tcs(tmp(2:end))-tcs(tmp(1:end-1)))./(tmptag(:,jch)+(tmptag(:,jch)==0))];
                    else
                        tcs = filter(tcs,jch:nch:end);
                        tcs = cumsum([zeros(1,size(tcs,2)); tcs]);
                        tim = [tim; (tcs(tmp(2:end),:)-tcs(tmp(1:end-1),:))];
                    end
                end
                in_line = false;

                tmpy = [rsty];
                tmpf = [rstf];
                tmpt = [rsttcspc];
                tmpm = [rstmarker];

                rsty      = [];
                rstf      = [];
                rsttcspc  = [];
                rstmarker = [];
            end

            tmp = find(tmpm~=0);
        end;
    end

    close(h);

    tag = reshape(tag, para(2), para(3), nch);
    tag(:, 2:2:end, :)  = tag(end:-1:1, 2:2:end, :);
    tag = permute(tag, [2 1 3]);

    nx        = head.ScanWidthX;
    ny        = head.ScanWidthY;

    if ~isempty(tim)&&(size(tim,1)==nx*ny)
        tim = reshape(tim, para(2), para(3), numel(tim)/prod(para(2:3))/nch, nch);
        tim(:, 2:2:end, :, :) = tim(end:-1:1, 2:2:end, :, :);
        tim = permute(tim, [2 1 3:ndims(tim)]);
        tim = squeeze(tim);
    end;


    figure;
    x = head.ScanStartX+(1:size(tag,1))*head.ScanResolution;
    y = head.ScanStartY+(1:size(tag,2))*head.ScanResolution;
    imagesc(y,x,sum(tag,3));
    set(gca,'DataAspectRatio', [1,1,1]);
    set(gca,'XDir','reverse');
    set(gca,'YDir','reverse');
    xlabel('x / µm');
    ylabel('y / µm');
end;
