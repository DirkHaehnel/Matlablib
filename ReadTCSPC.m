function [bin, tcspcdata, head] = ReadTCSPC(name)

% [bin, tcspcdata, head] = ReadTCSPC(name)
% name : file name
% (c) Christoph Pieper (2010)

if strcmp(name(end-2:end),'pt3')
    head = pt3v2Read(name);
    NChannels  = 2^12;
    NCounts   = head.Records;
elseif strcmp(name(end-2:end),'t3r')
    head = tttrRead(name);
    NChannels  = 2^12;
    NCounts   = head.Records;
elseif strcmp(name(end-2:end),'ht3')
	head = ht3v2read(name);
    if isempty(head)
        if nargin<3 || isempty(para)
            Timeunit = 50e-9;
            Resolution = 2e-12;
        else
            Timeunit = para(1)*1e-9;
            Resolution = para(2)*1e-9;
        end
        NChannels  = ceil(Timeunit/Resolution);
        tmp = dir(name);
        NCounts = tmp.bytes/4;
    else
        NChannels = 2^15;
        tmp = dir(name);
        NCounts = tmp.bytes/4;
    end
else
    disp('Not a valid file!');
    return
end

bin = [];
tcspcdata = [];
tcspcfile = [name(1:end-4) '.tcspc'];
matfile = [name(1:end-4) '.mat'];

if ( exist(matfile, 'file') == 2 )
    load(matfile);
    if ((~isempty(who('-file', matfile, 'UsedBunchs'))) && (isfield(res, 'tcspcdata2')))
        tcspcdata = sum(res.tcspcdata2(:,:,UsedBunchs),3);
    else
        tcspcdata = res.tcspcdata;
    end
elseif ( exist(tcspcfile, 'file') == 2) % check if tcspc data was saved
    load(tcspcfile, '-mat');
else
    photons = 1e6; % number of photons processed one at a time

    tst = 1;
    cnt = 0;
    bin = 0:NChannels-1;
    h = waitbar(0,'Please wait...');
    dind = [];
    rev_dind = zeros(8,1);
    while tst>0
        if strcmp(name(end-2:end),'pt3')
            [tmpy, tmpx, flv, bla, tst] = pt3v2Read(name, [cnt+1, photons]);
        elseif strcmp(name(end-2:end),'t3r')
            [tmpy, tmpx, flv, bla, tst] = tttrRead(name, [cnt+1, photons]);
        elseif strcmp(name(end-2:end),'ht3')
            if isempty(head)
                [tmpy, tmpx, flv, special, tst, overcount, head] = ht3read(name, [cnt+1, photons]);
            else
                [tmpy, tmpx, flv, special, tst, head] = ht3v2read(name, [cnt+1, photons]);
            end
        end
        if isempty(dind)
            %dind = union(flv,flv);
            dind = unique(flv);
            dnum = length(dind);
            tcspcdata = zeros(NChannels,dnum);
            for j=1:length(dind)
                rev_dind(dind(j)+1) = j;
            end
        end
        cnt = cnt+tst;
        waitbar(cnt/NCounts,h)
        if (tst>0) && (~isempty(tmpy))
            for j=1:length(flv)
                dd = rev_dind(flv(j)+1);
                if (dd > 0)
                    tcspcdata(tmpx(j)+1, dd) = tcspcdata(tmpx(j)+1, dd) + 1;
                end
            end
        end
    end
    close(h)

    ind = sum(tcspcdata,2)==0;
    bin(ind) = [];
    tcspcdata(ind,:) = [];
    bin = bin';
    
    save(tcspcfile, 'bin', 'tcspcdata');
end

end
