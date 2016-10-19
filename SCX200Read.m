function [tag, tim, tch, bin, head, tau, filter] = SCX200Read(name, filter)

[tmp, tmp, tmp, head] = tttrRead(name,[1 1]);

h = waitbar(0,'Please wait ...'); drawnow

para = [head.ScanTimePerPix/head.GlobClock/1e-3, head.ScanWidthX + head.ScanPause, head.ScanWidthY];

maxCounts = 2^18;
blocks = floor(head.NCounts/maxCounts)-1;

waitbar(1/max([1 blocks]));
[y, flag, tcspc, head, num, overcount, marker] = tttrRead(name, [1 head.NCounts - blocks*maxCounts]);
y(marker>0) = [];
flag(marker>0) = [];
tcspc(marker>0) = [];        

cnt = head.NCounts - blocks*maxCounts;    
my = y(1);

tmp = y>=floor(max(y)/para(1))*para(1);
rsty = y(tmp);
rstf = flag(tmp);
rsttcspc = tcspc(tmp);
tmpy = y(1:end-length(rsty));
tmpf = flag(1:end-length(rsty));
tmpt = tcspc(1:end-length(rsty));
tmptag = tttr2bin(tmpy,para(1),tmpf);
nch = size(tmptag,2);
tag = tmptag;
bin = min(tcspc):max(tcspc);

tau = [];
if nargin>1 & ~isempty(filter) & ~(size(filter,1)==length(bin))
    tau = filter;
    maxtau = max(tau);
    if length(tau)>1
        filter = exp(-head.Resolution*bin'*(1./tau(:)'));
        tst = 1; cc = 0; 
        weight = zeros(length(bin),nch);
        while tst>0
            [tmp, flv, tmpx, tmp, tst, tmp, marker] = tttrRead(name,[cc+1 maxCounts]);        
            tmpy(marker>0) = [];
            flv(marker>0) = [];
            tmpx(marker>0) = [];        
            if tst>0
                cc = cc+tst;
                for jch=1:nch
                    weight(:,jch) = weight(:,jch) + mhist(tmpx(flv==jch-1), bin);
                end
            end
        end
        filter = reshape(repmat(filter,[nch 1]),length(bin),size(filter,2)*nch);
        for jch=1:nch
            [m,m] = max(weight(:,jch));
            filter(1:m+9,jch:nch:end) = 0;
            weight(1:m+9,jch) = 0;        
            m = max(bin(filter(:,jch)>0));
            filter(m-9:end,jch:nch:end) = 0;
            weight(m-9:end,jch) = 0;        
        end
            filter = filter./(ones(size(filter,1),1)*sum(filter));
        for jch = 1:nch
            ind = weight(:,jch)>0;
            filter(ind,jch:nch:end) = ((filter(ind,jch:nch:end)'*diag(1./weight(ind,jch))*filter(ind,jch:nch:end))\...
                (diag(1./weight(ind,jch))*filter(ind,jch:nch:end))')';
        end
    else
        filter = zeros(2,nch);
        for jch=1:nch
            [m,m] = max(mhist(tmpt(tmpf==jch-1), bin));
            tmp = round(3*maxtau/head.Resolution)+m+9;
            filter(:,jch) = [bin(m+9); bin(min([end-10 tmp]))];
        end
    end
end

if nargin==1
    for jch=1:nch
        [m,m] = max(mhist(tmpt(tmpf==jch-1), bin));
        filter(:,jch) = [bin(m+9); bin(end-10)];
    end
end

for jch=1:nch
    tcs = tmpt(tmpf==jch-1);
    if nargout>3
        tch(:,jch) = mhist(tcs,bin);
    end
    tmp = cumsum([1; tag(:,jch)]);
    if size(filter,1)==2
        tst = tcs>filter(1,jch) & tcs<=filter(2,jch);
        tcs = cumsum([0; (tcs - filter(1,jch)).*tst]);
        tcs = tcs(tmp(2:end))-tcs(tmp(1:end-1));
        tst = cumsum([0; tst]);
        tst = tst(tmp(2:end))-tst(tmp(1:end-1));
        tim(:,jch) = tcs./(tst+(tst==0));
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
    y(marker>0) = [];
    flag(marker>0) = [];
    tcspc(marker>0) = [];
    cnt = cnt + num;
    y = [rsty; overcount + y];
    overcount = overcount + tmp;
    flag = [rstf; flag];
    tcspc = [rsttcspc; tcspc];
    tmp = y>=floor(max(y)/para(1))*para(1);
    rsty = y(tmp);
    rstf = flag(tmp);
    rsttcspc = tcspc(tmp);
    tmpy = y(1:end-length(rsty));
    tmpf = flag(1:end-length(rsty));
    tmpt = tcspc(1:end-length(rsty));
    ind = size(tag,1);
    tmptag = tttr2bin([my+para(1)*ind; tmpy],para(1),[0; tmpf]);
    tmptag(1,1) = tmptag(1,1) - 1;
    tag(ind+1:ind+size(tmptag,1),:) = tmptag;
    for jch=1:nch
        tcs = tmpt(tmpf==jch-1);
        if nargout>3
            tch(:,jch) = tch(:,jch) + mhist(tcs,bin);
        end
        tmp = cumsum([1; tmptag(:,jch)]);
        if size(filter,1)==2
            tst = tcs>filter(1,jch) & tcs<=filter(2,jch);
            tcs = cumsum([0; (tcs - filter(1,jch)).*tst]);
            tcs = tcs(tmp(2:end))-tcs(tmp(1:end-1));
            tst = cumsum([0; tst]);
            tst = tst(tmp(2:end))-tst(tmp(1:end-1));
            tim(ind+1:ind+size(tmptag,1),jch) = tcs./(tst+(tst==0));
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

t = reshape(1:para(2)*para(3),para(2),para(3));
t(:,2:2:end) = t(end:-1:1,2:2:end);
t = reshape(t,1,para(2)*para(3));

tmptag = sum(tag,2);
j = 0;
tmp = tmptag(j+t(1:para(2)*(para(3)-1)))'*tmptag(j+t((para(2)+1):para(2)*para(3)));
for j=5:5:min([50, length(tmptag)-max(t)])
    tmp = [tmp, tmptag(j+t(1:para(2)*(para(3)-1)))'*tmptag(j+t((para(2)+1):para(2)*para(3)))]; 
end
[pos, pos] = max(tmp);      
tmp = tmp(max([1,pos-1]));
for j=5*max([0,pos-2])+1:5*max([0,pos-2])+10
    tmp = [tmp, tmptag(j+t(1:para(2)*(para(3)-1)))'*tmptag(j+t((para(2)+1):para(2)*para(3)))]; 
end
bini = 5*max([0,pos-2]):0.1:5*max([0,pos-2])+10;
tmpi = interp1(5*max([0,pos-2]):5*max([0,pos-2])+10,tmp,bini,'cubic');
[tmpi, tmpi] = max(tmpi);      
offset = bini(tmpi);

if offset==0
    para(3) = para(3) - 1;
    head.ScanWidthY = para(3);
    if nch>1
        tag(1:para(2)-25,:) = [];
        tim(1:para(2)-25,:,:) = [];    
        tmptag = sum(tag')';
    else
        tag(1:para(2)-25) = [];
        tim(1:para(2)-25,:) = [];    
        tmptag = tag;
    end
    j = 0;
    tmp = tmptag(j+t(1:para(2)*(para(3)-1)))'*tmptag(j+t((para(2)+1):para(2)*para(3)));
    for j=5:5:50 
        tmp = [tmp, tmptag(j+t(1:para(2)*(para(3)-1)))'*tmptag(j+t((para(2)+1):para(2)*para(3)))]; 
    end
    [pos, pos] = max(tmp);      
    tmp = tmp(max([1,pos-1]));
    for j=5*max([0,pos-2])+1:5*max([0,pos-2])+10
        tmp = [tmp, tmptag(j+t(1:para(2)*(para(3)-1)))'*tmptag(j+t((para(2)+1):para(2)*para(3)))]; 
    end
    bini = 5*max([0,pos-2]):0.1:5*max([0,pos-2])+10;
    tmpi = interp1(5*max([0,pos-2]):5*max([0,pos-2])+10,tmp,bini,'cubic');
    [tmpi, tmpi] = max(tmpi);      
    offset = bini(tmpi);
    head = setfield(head, 'FigureOffset', offset-25);    
else
    head = setfield(head, 'FigureOffset', offset);
end

tag = (floor(offset)+1-offset)*tag(floor(offset)+(1:para(2)*para(3)),:) + (offset-floor(offset))*tag(ceil(offset)+(1:para(2)*para(3)),:);
tim = (floor(offset)+1-offset)*tim(floor(offset)+(1:para(2)*para(3)),:,:) + (offset-floor(offset))*tim(ceil(offset)+(1:para(2)*para(3)),:,:);
tag = reshape(tag,para(2),para(3),nch);
tag(:,2:2:end,:) = tag(end:-1:1,2:2:end,:);
tim = reshape(tim,para(2),para(3),prod(size(tim))/prod(para(2:3))/nch,nch);   
tim(:,2:2:end,:,:) = tim(end:-1:1,2:2:end,:,:);   
tim = squeeze(tim);
if size(filter,1)==2
    tim = tim*head.Resolution;
end
close(h);

