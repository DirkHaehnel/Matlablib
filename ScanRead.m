function [tag, tim, tcspc, head, thist, v, raw] = ScanRead(name, para)

% The function [tag, tim, tcspc, head, thist, v, raw] = ScanRead(name, para) reads image data from file 'name' 
% para(1) = bin time in 100 nsec
% para(2:3) = image size
% tag = intensity image
% tim = average lifetime image
% tcspc = photon lifetime data
% head = measurement parameters
% raw = raw image, without alignment

peakshift = 15;
fin = fopen(name,'r');
if (fin==-1)
    errordlg('Cannot open specified file. Please try again.');
else        
    h = waitbar(0,'Please wait ...'); drawnow
    fseek(fin, 48+256, 0);
    head = struct('NChannels', fread(fin, 1, 'int32'));
    head = setfield(head, 'NCurves', fread(fin, 1, 'int32'));
    head = setfield(head, 'ActiveCurve', fread(fin, 1, 'int32'));
    head = setfield(head, 'MeasMode', fread(fin, 1, 'int32'));
    head = setfield(head, 'HistMode', fread(fin, 1, 'int32'));
    head = setfield(head, 'Range', fread(fin, 1, 'int32'));
    head = setfield(head, 'Offset', fread(fin, 1, 'int32'));
    head = setfield(head, 'AcquTime', fread(fin, 1, 'int32'));
    head = setfield(head, 'LinLog', fread(fin, 1, 'int32'));
    head = setfield(head, 'MinAx', fread(fin, 1, 'int32'));
    head = setfield(head, 'MaxAx', fread(fin, 1, 'int32'));
    head = setfield(head, 'MinAxCnt', fread(fin, 1, 'int32'));
    head = setfield(head, 'MaxAxCnt', fread(fin, 1, 'int32'));
    head = setfield(head, 'CFD0', fread(fin, 1, 'int32'));
    head = setfield(head, 'CFDmin', fread(fin, 1, 'int32'));
    head = setfield(head, 'Sync', fread(fin, 1, 'int32'));
    head = setfield(head, 'Resolution', fread(fin, 1, 'float'));
    head = setfield(head, 'GlobClock', fread(fin, 1, 'int32'));
    head = setfield(head, 'NCounts', fread(fin, 1, 'uint32'));
    
    waitbar(0.25);
    maxCounts = 2^20;
    blocks = floor(head.NCounts/maxCounts);
    
    y = fread(fin, head.NCounts - blocks*maxCounts, 'uint32');   
    flag = floor(y/2^30);
    y = y - floor(y/2^28)*2^28;
    tcspc = floor(y/2^16);
    y = y - tcspc*2^16;
    tcspc = 2^12 + 1 - tcspc; 
    y = y + 2^16*cumsum(~(flag==1));
    cnt = sum(~(flag==1));
    y(~(flag==1))=[];
    tcspc(~(flag==1))=[];
    
    y = [y(1); y];
    tcspc = [max(tcspc); tcspc];
    
    rst = y(y>floor(max(y)/para(1))*para(1));
    rsttcspc = tcspc(end-length(rst)+1:end);
    y = y(1:end-length(rst));
    tcspc = tcspc(1:end-length(rst));
    
    tag = tttr2bin(y,para(1));
    v = min(tcspc):max(tcspc);
    thist = mhist(tcspc,v);
    [peakpos, peakpos] = max(thist);
    v(peakpos);
    tim = tcspc;
    tmpy = y;
    tmpy(tim<v(peakpos+peakshift))=[];
    tim(tim<v(peakpos+peakshift))=[];
    tmpy = [y(1); tmpy];
    tim = [tcspc(1); tim];
    tmpy = tttr2bin(tmpy,para(1));
    tmp = cumsum([1; tmpy]);
    tim = cumsum([0; tim-v(peakpos)]);
    tim = (tim(tmp(2:end))-tim(tmp(1:end-1)))./(tmpy+(tmpy==0));
    
    for j=1:blocks
        y = fread(fin, maxCounts, 'uint32');   
        flag = floor(y/2^30);
        y = y - floor(y/2^28)*2^28;
        tcspc = floor(y/2^16);
        y = y - tcspc*2^16;
        tcspc = 2^12 + 1 - tcspc; 
        y = y + 2^16*(cnt+cumsum(~(flag==1)));
        cnt = cnt + sum(~(flag==1));
        y(~(flag==1))=[];
        tcspc(~(flag==1))=[];
        
        y = [rst; y];
        tcspc = [rsttcspc; tcspc];
        rst = y(y>floor(max(y)/para(1))*para(1));
        rsttcspc = tcspc(end-length(rst)+1:end);
        y = y(1:end-length(rst));
        tcspc = tcspc(1:end-length(rst));
        
        thist = thist + mhist(tcspc,v);
        tmptim = tcspc;
        tmpy = y;
        tmpy(tmptim<v(peakpos+peakshift))=[];
        tmptim(tmptim<v(peakpos+peakshift))=[];
        tmpy = tttr2bin(tmpy,para(1));
        tmp = cumsum([1; tmpy]);
        tmptim = cumsum([0; tmptim-v(peakpos)]);
        tim = [tim; (tmptim(tmp(2:end))-tmptim(tmp(1:end-1)))./(tmpy+(tmpy==0))];
        tag = [tag; tttr2bin(y,para(1))];
    end
    
    fclose(fin);
    
    waitbar(0.5);
    t = reshape(1:para(2)*para(3),para(2),para(3));
    t(:,2:2:end) = t(end:-1:1,2:2:end);
    t = reshape(t,1,para(2)*para(3));
    j = 0;
    tmp = tag(j+t(1:para(2)*(para(3)-1)))'*tag(j+t((para(2)+1):para(2)*para(3)));
    for j=1:20
        tmp = [tmp, tag(j+t(1:para(2)*(para(3)-1)))'*tag(j+t((para(2)+1):para(2)*para(3)))]; 
    end
    [offest, offset] = max(tmp);
    offset = offset - 1;
    raw = reshape(tag(1:para(2)*para(3)),para(2),para(3));
    raw(:,2:2:end) = raw(end:-1:1,2:2:end);
    tcspc(1:sum(tag(1:offset))) = [];
    tag(1:offset) = [];
    waitbar(0.75);
    tag = reshape(tag(1:para(2)*para(3)),para(2),para(3));
    tag(:,2:2:end) = tag(end:-1:1,2:2:end);
    tim(1:offset) = [];
    tim = reshape(tim(1:para(2)*para(3)),para(2),para(3));   
    tim(:,2:2:end) = tim(end:-1:1,2:2:end);   
    tim = tim*head.Resolution;
    close(h);
    if nargin==3 & ~isempty(flag)
        subplot(121); imagesc(log10(tag)); axis off; colorbar('h'); axis image; 
        subplot(122); imagesc(tim); axis image; axis off; colorbar('h'); 
        colormap hot
    end
end

% subplot(121); imagesc(log10(tag)); colorbar('h'); axis image; subplot(122); imagesc(tim); axis image; colorbar('h'); colormap hot

