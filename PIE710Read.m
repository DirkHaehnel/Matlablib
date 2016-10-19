function [tag, tim, tch, bin, head] = PIE710Read(name, filter)

[tmp, tmp, tmp, head] = tttrRead(name,[1 1]);

h = waitbar(0,'Please wait ...'); drawnow

para = [head.ScanTimePerPix*2e3, head.ScanWidthX, head.ScanWidthY];

maxCounts = 2^19;
blocks = floor(head.NCounts/maxCounts)-1;

waitbar(1/blocks);
[y, flag, tcspc, head, num, overcount, marker] = tttrRead(name, [1 head.NCounts - blocks*maxCounts]);
marker = bitget(marker,3);

cnt = head.NCounts - blocks*maxCounts;    
my = y(1);

tmp = y>=floor(max(y)/para(1))*para(1);
rsty = y(tmp);
rstf = flag(tmp);
rsttcspc = tcspc(tmp);
rstmarker = marker(tmp);
tmpy = y(1:end-length(rsty));
tmpf = flag(1:end-length(rsty));
tmpt = tcspc(1:end-length(rsty));
tmpm = marker(1:end-length(rsty));
tag = tttr2bin(tmpy,para(1),tmpf);
blip = tttr2bin(tmpy,para(1),tmpm==1);
nch = size(tag,2);
bin = min(tcspc):max(tcspc);

for jch=1:nch
    tcs = tmpt(tmpf==jch-1);
    if nargout>2
        tch(:,jch) = mhist(tcs,bin);
    end
    tmp = cumsum([1; tag(:,jch)]);
    if nargin<2
        tcs = cumsum([0; tcs]);
        tim(:,jch) = (tcs(tmp(2:end))-tcs(tmp(1:end-1)))./(tag(:,jch)+(tag(:,jch)==0));
    else
        tcs(tcs<bin(1)) = bin(1);
        tcs(tcs>bin(end)) = bin(end);
        tcs = filter(tcs-bin(1)+1,jch:nch:end);
        tst = abs(tcs)>0;
        tcs = cumsum([zeros(1,size(tcs,2)); tcs]);
        tcs = tcs(tmp(2:end),:)-tcs(tmp(1:end-1),:);
        tst = cumsum([zeros(1,size(tcs,2)); tst]);
        tst = tst(tmp(2:end),:)-tst(tmp(1:end-1),:);
        tim(:,:,jch) = tcs;%./(tst+(tst==0));
    end            
end

for j=1:blocks
    waitbar(j/blocks);
    [y, flag, tcspc, head, num, tmp, marker] = tttrRead(name, [cnt+1 maxCounts]);
    marker = bitget(marker,3);
    cnt = cnt + num;
    y = [rsty; overcount + y];
    overcount = overcount + tmp;
    flag = [rstf; flag];
    tcspc = [rsttcspc; tcspc];
    marker = [rstmarker; marker];
    tmp = y>=floor(max(y)/para(1))*para(1);
    rsty = y(tmp);   
    rstf = flag(tmp);   
    rsttcspc = tcspc(tmp);        
    rstmarker = marker(tmp);
    tmpy = y(1:end-length(rsty));
    tmpf = flag(1:end-length(rsty));
    tmpt = tcspc(1:end-length(rsty));
    tmpm = marker(1:end-length(rsty));
    ind = size(tag,1);
    tmptag = tttr2bin([my+para(1)*ind; tmpy],para(1),[0; tmpf]);    
    tmptag(1,1) = tmptag(1,1) - 1;
    tag(ind+1:ind+size(tmptag,1),:) = tmptag;
    tmpm = tttr2bin([my+para(1)*ind; tmpy],para(1),[true; tmpm==1]);    
    tmpm(1) = tmpm(1) - 1;
    blip(ind+1:ind+size(tmptag,1),:) = tmpm;
    for jch=1:nch
        tcs = tmpt(tmpf==jch-1);
        if nargout>2
            tch(:,jch) = tch(:,jch) + mhist(tcs,bin);
        end
        tmp = cumsum([1; tmptag(:,jch)]);
        if nargin<2
            tcs = cumsum([0; tcs]);
            tmp = cumsum([1; tmptag(:,jch)]);
            tim(ind+1:ind+size(tmptag,1),jch) = (tcs(tmp(2:end))-tcs(tmp(1:end-1)))./(tmptag(:,jch)+(tmptag(:,jch)==0));
        else
            tcs(tcs<bin(1)) = bin(1);
            tcs(tcs>bin(end)) = bin(end);
            tcs = filter(tcs-bin(1)+1,jch:nch:end);
            tst = abs(tcs)>0;
            tcs = cumsum([zeros(1,size(tcs,2)); tcs]);
            tcs = tcs(tmp(2:end),:)-tcs(tmp(1:end-1),:);
            tst = cumsum([zeros(1,size(tcs,2)); tst]);
            tst = tst(tmp(2:end),:)-tst(tmp(1:end-1),:);
            tim(ind+1:ind+size(tmptag,1),:,jch) = tcs;%./(tst+(tst==0));
        end
    end    
end

tmptag = tag; tmptim = tim;
t = 1:size(tmptag,1);
t = t(blip==1);
dt = mean(diff(t));
tag = []; tim = [];
for j=1:length(t)
    tag(:,j,:) = tmptag(round(t(j)+dt*head.ScanTStartTo):round(t(j)+dt*head.ScanTStopTo),:);
    tim(:,j,:,:) = tmptim(round(t(j)+dt*head.ScanTStartTo):round(t(j)+dt*head.ScanTStopTo),:,:);    
end
tag(:,2:2:end,:,:) = tag(end:-1:1,2:2:end,:,:);
tim(:,2:2:end,:,:) = tim(end:-1:1,2:2:end,:,:);
tim = squeeze(tim);
close(h);
