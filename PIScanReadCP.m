function [intens, flim, tcspc, head] = PIScanReadCP(name, flimtype, readFrames, LSMReductionFactor)

% Parameters
% name :        Filename
% flimtype :    0 = default, 1 = no flim / only intens, 2 = histogram
% readFrames:   number of frames to read (LSM files only)
% LSMRedFactor: reduce TCSPC bins by 2^x

head = ht3v2read_head(name);
photons = 1e6; % number of photons processed one at a time

if nargin<4
    LSMReductionFactor = 2;
    if (strcmp(head.Ident,'HydraHarp'))
        LSMReductionFactor = 5;
    end
end
if nargin<3
    readFrames = inf;
end
if nargin<2
    flimtype = 0;
end

isMTSImage = (regexp(head.CommentField, 'MT Scan Image') == 1);
isLSMImage = 0;
if (~isempty(head.ImgHdr))
    isLSMImage = (head.ImgHdr.Ident == 3);
else
    if (isMTSImage)
        
        % Picture Size
        ps = regexp(head.CommentField, 'Picture Size: [0-9]+\.[0-9]+ um x [0-9]+\.[0-9]+', 'match');
        sizes = regexp(ps, '(?<SizeX>[0-9]+\.[0-9]+) um x (?<SizeY>[0-9]+\.[0-9]+)', 'tokens');
        head.ImgHdr.SizeX = str2num(sizes{1}{1}{1});
        head.ImgHdr.SizeY = str2num(sizes{1}{1}{2});
        
        %Picture Anz Pixel
        ps = regexp(head.CommentField, 'Pixels: [0-9]+ x [0-9]+', 'match');
        sizes = regexp(ps, '(?<PixelX>[0-9]+) x (?<PixelY>[0-9]+)', 'tokens');
        head.ImgHdr.PixX = str2num(sizes{1}{1}{1});
        head.ImgHdr.PixY = str2num(sizes{1}{1}{2});

        %Picture UL Corner
        ps = regexp(head.CommentField, 'Upper Left \(x:y:z\) \[um\]: [0-9]+\.[0-9]+ : [0-9]+\.[0-9]+ : [0-9]+\.[0-9]+', 'match');
        sizes = regexp(ps, '(?<SizeX>[0-9]+\.[0-9]+) : (?<SizeY>[0-9]+\.[0-9]+) : (?<SizeY>[0-9]+\.[0-9]+)', 'tokens');
        head.ImgHdr.X0 = str2num(sizes{1}{1}{1});
        head.ImgHdr.Y0 = str2num(sizes{1}{1}{2});
    else
        intens = [];
        flim = [];
        return
    end
end

[sync, tcspc, chan, special, num, overcount, head2] = ht3v2read_old(name, [1 1e5]);

dind = unique(chan(special == 0));
ndet = length(dind);

rev_dind = zeros(8,1);
for j=1:length(dind)
    rev_dind(dind(j)+1) = j;
end

if (isMTSImage)
    intens = zeros(head.ImgHdr.PixY, head.ImgHdr.PixX, ndet);
else
    intens = zeros(head.ImgHdr.PixY, head.ImgHdr.PixX);    
end

flim = [];
if (isLSMImage)
    flim = zeros(head.ImgHdr.PixY,head.ImgHdr.PixX,2^(12-LSMReductionFactor));
else
    if (flimtype == 2)
        if (strcmp(head.Ident,'HydraHarp'))
            flim = zeros(head.ImgHdr.PixY,head.ImgHdr.PixX,2^(15-LSMReductionFactor));
        else
            flim = zeros(head.ImgHdr.PixY,head.ImgHdr.PixX,2^(12-LSMReductionFactor));
        end   
    else
        flim = cell(head.ImgHdr.PixY,head.ImgHdr.PixX);
    end
end

cnt = 1;
num = photons;

h = waitbar(0,'Reading image data');

if (isMTSImage)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % monodirectional (MTScan LabView)
    Turns = [];
    dT = Inf;

    currX = 0;
    currY = 1;
    tStart = 0;
    Pixels = head.ImgHdr.PixX * head.ImgHdr.PixY;
    
    tcspc2 = [];

    while (num == photons)
        [sync, tcspc, chan, special, num, overcount, head2] = ht3v2read_old(name, [cnt photons]);
        cnt = cnt + num;
        
        if (length(Turns) < 2)
            Turns2 = sync(chan==4 & special==1);
            Turns = [Turns; Turns2];
        end
        if (length(Turns) >= 2)
            dT = Turns(2) - Turns(1);
        end
        trigposs = find(chan==4 & special==1);

        % photons that belong to last pixel
        if (currX + currY > 1)
            t2 = tStart + dT;
            
            if (length(trigposs) >= 1)
                tmptcspc = tcspc(1:trigposs(1)-1);
                tmpy = sync(1:trigposs(1)-1);
                tmpchan = chan(1:trigposs(1)-1);
            else
                tmptcspc = tcspc;
                tmpy = sync;
                tmpchan = chan;                    
            end

            if (currX == head.ImgHdr.PixX)
                idx1 = find(tmpy>=tStart , 1, 'first');
                idx2 = find(tmpy<=t2 , 1, 'last');
            else
                idx1 = 1;
                if (length(trigposs) >= 1)
                    idx2 = trigposs(1)-1;
                else
                    idx2 = length(tmptcspc);
                end
            end
            tmptcspc = tmptcspc(idx1:idx2);
            tmpchan = tmpchan(idx1:idx2);          

            if (flimtype == 1)
                for di = 1:ndet 
                    intens(currY, currX, di) = intens(currY, currX, di) + length(tmptcspc(tmpchan==dind(di)));
                end
            elseif (flimtype == 2)
                for k=1:length(tmptcspc)
                    channl = bitshift(tmptcspc(k),-LSMReductionFactor) + 1;
                    flim(currY, currX, channl) = flim(currY, currX, channl) + 1;
                    intens(currY, currX, rev_dind(tmpchan(k))) = intens(currY, currX, rev_dind(tmpchan(k))) + 1;
                end
            else
                flim{currY, currX} = [flim{currY, currX}; tmptcspc];
                for di = 1:ndet 
                    intens(currY, currX, di) = intens(currY, currX, di) + length(tmptcspc(tmpchan==dind(di)));
                end
            end
        end
        
        % photons for new pixels
        for j=1:length(trigposs)
            currX = currX + 1;
            if (currX > head.ImgHdr.PixX)
                currX = 1;
                currY = currY + 1;
            end
            tStart = sync(trigposs(j));
            
            if (j == length(trigposs))
                tmpy = sync(trigposs(j)+1:end);
            else
                tmpy = sync(trigposs(j)+1:trigposs(j+1)-1);
            end

            idx1 = trigposs(j)+1;
            if (j == length(trigposs)) || (currX == head.ImgHdr.PixX)
                t2 = tStart + dT;
                idx2 = trigposs(j) + find(tmpy<=t2 , 1, 'last');
            else
                idx2 = trigposs(j+1)-1;
            end
            
            tmptcspc = tcspc(idx1:idx2);
            tmpchan = chan(idx1:idx2);
                         
            if (flimtype == 1)
                for di = 1:ndet 
                    intens(currY, currX, di) = intens(currY, currX, di) + length(tmptcspc(tmpchan==dind(di)));
                end
            elseif (flimtype == 2)
                for k=1:length(tmptcspc)
                    channl = bitshift(tmptcspc(k),-LSMReductionFactor) + 1;
                    flim(currY, currX, channl) = flim(currY, currX, channl) + 1;
                    intens(currY, currX, rev_dind(tmpchan(k))) = intens(currY, currX, rev_dind(tmpchan(k))) + 1;
                end
            else
                flim{currY, currX} = [flim{currY, currX}; tmptcspc];
                for di = 1:ndet 
                    intens(currY, currX, di) = intens(currY, currX, di) + length(tmptcspc(tmpchan==dind(di)));
                end
            end            

            waitbar((head.ImgHdr.PixX*(currY-1)+currX)/Pixels,h); drawnow;        
        end

        tcspc(special == 1) = [];
        chan(special == 1) = [];
        if (strcmp(head.Ident,'HydraHarp'))
            if (isempty(tcspc2))
                tcspc2 = zeros(ndet2^15,ndet);
            end
            for di = 1:ndet
                tcspc2(:,di) = tcspc2(:,di) + mHist(tcspc(chan == dind(di)),0:2^15-1);
            end
        else
            if (isempty(tcspc2))
                tcspc2 = zeros(2^12,ndet);
            end
            for di = 1:ndet
                tcspc2(:,di) = tcspc2(:,di) + mHist(tcspc(chan == dind(di)),0:2^12-1);
            end
        end    
    end
    tcspc = tcspc2;
elseif (isLSMImage)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % LSM image monodirectional
    Turns = [];
    dT = Inf;

    currX = 0;
    currY = 1;
    currFrame = 1;
    tStart = 0;
    Pixels = head.ImgHdr.PixX * head.ImgHdr.PixY;
    trigLStart = head.ImgHdr.LineStart;
    trigLStop = head.ImgHdr.LineStop;
    trigFrame = head.ImgHdr.Frame;
    
    tcspc2 = [];
    sync = [];
    tcspc = [];
    chan = [];
    special = [];
    addfrom = 0;
    overcount = 0;

    while ((num == photons) && (currFrame <= readFrames))
        [sync_tmp, tcspc_tmp, chan_tmp, special_tmp, num, overcount2, head2] = ht3v2read_old(name, [cnt photons]);
        if (addfrom > 0)
            sync = [sync(addfrom:end); sync_tmp + overcount];
            tcspc = [tcspc(addfrom:end); tcspc_tmp];
            chan = [chan(addfrom:end); chan_tmp];
            special = [special(addfrom:end); special_tmp];
        else
            sync = sync_tmp + overcount;
            tcspc = tcspc_tmp;
            chan = chan_tmp;
            special = special_tmp;
        end
        overcount = overcount2 + overcount;
        cnt = cnt + num;
        waitbar(cnt/head.Records,h); drawnow;

        trigpossLStart = find(chan==trigLStart & special==1);
        trigpossLStop = find(chan==trigLStop & special==1);
        trigpossFrame = find(chan==4 & special==1);
        if (currFrame == 1)
            trigpossFrame = [0; trigpossFrame];
        end
        if (addfrom == 0)
            trigpossLStop = trigpossLStop(2:end);
        end
        
        numFrames = length(trigpossFrame)-1;
        
        for i=1:numFrames
            trigpossLStart2 = trigpossLStart(trigpossLStart > trigpossFrame(i) & trigpossLStart < trigpossFrame(i+1));
            trigpossLStop2 = trigpossLStop(trigpossLStop > trigpossFrame(i) & trigpossLStop < trigpossFrame(i+1));
            numPix = length(trigpossLStop2);
            
            % Splitup into Lines
            for j=1:numPix
                dT = (sync(trigpossLStop2(j)) - sync(trigpossLStart2(j))) / head.ImgHdr.PixX;
                t1 = sync(trigpossLStart2(j));
                t2 = t1 + dT;
                
                tt_tcspc = tcspc(trigpossLStart2(j)+1:trigpossLStop2(j)-1);
                tt_sync = sync(trigpossLStart2(j)+1:trigpossLStop2(j)-1);
                tt_photons = length(tt_sync);
                % Splitup into Pixel
                currX = 1;
                for k=1:tt_photons
                    while (t2 < tt_sync(k))
                        t2 = t2 + dT;
                        currX = currX  + 1;
                    end
                    if (flimtype ~= 1)
                        channl = bitshift(tt_tcspc(k),-LSMReductionFactor) + 1;
                        flim(j, currX, channl) = flim(j, currX, channl) + 1;
                    end
                    intens(j, currX) = intens(j, currX) + 1;
                end
            end
            currFrame = currFrame + 1;
            if (currFrame > readFrames)
                break;
            end
        end
        
        if (numFrames > 0)
            addfrom = trigpossFrame(numFrames)+1;
        else
            addfrom = 1;
        end

        tcspc_tmp(special_tmp == 1) = [];
        if (strcmp(head.Ident,'HydraHarp'))
            if (isempty(tcspc2))
                tcspc2 = mHist(tcspc_tmp,0:2^15-1);
            else
                tcspc2 = tcspc2 + mHist(tcspc_tmp,0:2^15-1);
            end
        else
            if (isempty(tcspc2))
                tcspc2 = mHist(tcspc_tmp,0:2^12-1);
            else
                tcspc2 = tcspc2 + mHist(tcspc_tmp,0:2^12-1);
            end
        end    
    end
    tcspc = tcspc2;
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Symphotime...
    [sync, tcspc, chan, special, num, overcount, head] = ht3v2read_old(name, [1 inf]);
    Turns = sync(chan==4 & special==1);
    if head.ImgHdr.Pattern==0
        % monodirectional (SymPhoTime)
        for j=1:length(Turns)-1
            dT = (Turns(j+1)-Turns(j)); 
            t1 = Turns(j)+head.ImgHdr.TStartTo*dT;
            t2 = Turns(j)+head.ImgHdr.TStopTo*dT;
            tmp = sync(sync>=t1 & sync<=t2) - t1;
            intens(j,:) = tttr2bin([0; tmp; t2-t1-1], (t2-t1)/(head.ImgHdr.PixX));
            intens(j,[1 end]) = intens(j,[1 end]) - 1;
            tmp = tcspc(sync>=t1 & sync<=t2);
            cnt = cumsum([1 intens(j,:)]);
            for k=1:head.ImgHdr.PixX
                
                if cnt(k+1)>cnt(k)
                    if (flimtype == 1)
                        flim{j+1,k} = 0;
                    elseif (flimtype == 2)
                        for m=cnt(k):cnt(k+1)-1
                            channl = bitshift(tmp(m),-LSMReductionFactor) + 1;
                            flim(j+1, k, channl) = flim(j+1, k, channl) + 1;
                            intens(j+1, k) = intens(j+1, k) + 1;
                        end
                    else
                        flim{j+1,head.ImgHdr.PixX-k+1} = tmp(cnt(k):cnt(k+1)-1);
                    end
                else
                    if (flimtype ~= 2)
                        flim{j+1,k} = 0;
                    end
                end
            end
            intens(j+1,:) = intens(j+1,end);
            waitbar(j/length(Turns),h); drawnow;
        end
    else
        %bidirectional (SymPhoTime)
        for j=1:2:length(Turns)-1
            dT = (Turns(j+1)-Turns(j));
            t1 = Turns(j)+head.ImgHdr.TStartTo*dT;
            t2 = Turns(j)+head.ImgHdr.TStopTo*dT;
            tmp = sync(sync>=t1 & sync<=t2) - t1;
            intens(j,:) = tttr2bin([0; tmp; t2-t1-1], (t2-t1)/(head.ImgHdr.PixX));
            intens(j,[1 end]) = intens(j,[1 end]) - 1;
            tmp = tcspc(sync>=t1 & sync<=t2);
            cnt = cumsum([1 intens(j,:)]);
            for k=1:head.ImgHdr.PixX
                if cnt(k+1)>cnt(k)
                    flim{j,k} = tmp(cnt(k):cnt(k+1)-1);
                else
                    flim{j,k} = 0;
                end
            end
            t1 = Turns(j)+head.ImgHdr.TStartFro*dT;
            t2 = Turns(j)+head.ImgHdr.TStopFro*dT;
            tmp = sync(sync>=t1 & sync<=t2) - t1;
            intens(j+1,:) = tttr2bin([0; tmp; t2-t1-1], (t2-t1)/(head.ImgHdr.PixX));
            intens(j+1,[1 end]) = intens(j+1,[1 end]) - 1;
            tmp = tcspc(sync>=t1 & sync<=t2);
            cnt = cumsum([1 intens(j+1,:)]);
            for k=1:head.ImgHdr.PixX
                if cnt(k+1)>cnt(k)
                    flim{j+1,head.ImgHdr.PixX-k+1} = tmp(cnt(k):cnt(k+1)-1);
                else
                    flim{j+1,k} = 0;
                end
            end
            intens(j+1,:) = intens(j+1,end:-1:1);
            waitbar(j/length(Turns),h); drawnow;
        end
    end
    
    tcspc(special == 1) = [];
    if (strcmp(head.Ident,'HydraHarp'))
        tcspc = mHist(tcspc,0:2^15-1);
    else
        tcspc = mHist(tcspc,0:2^12-1);
    end   
end
last = find(sum(tcspc,2)>0,1,'last');
tcspc = tcspc(1:last,:);

close(h);


