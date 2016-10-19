function [bin, tcspcdata, head]=ht3TCSPC(name)

% This program gives the TCSPC Histogram for .ht3 files. name = filename
% This program needs HT3_Read.m
% (c) Narain Karedla (2013)

photons=1e6; % 1 million time at a time for processing.
if strcmp(name(end-2:end),'ht3')
    head=HT3_Read(name);
    Timeunit=1/(head.SyncRate);
    Resolution=head.Resolution*1e-9;
    NChannels=ceil(Timeunit/Resolution);
    NCounts=head.Records;
    
    
else
    disp('Not a valid file!')
    return
end
 h = waitbar(0,'Please wait...');
    bin = [];
    tcspcdata = [];
    tcspcfile = [name(1:end-7) '.ht3tcspc'];
    %matfile = [name(1:end-4) '.mat'];

%     if ( exist(matfile, 'file') == 2 )
%         load(matfile);
%     if ((~isempty(who('-file', matfile, 'UsedBunchs'))) && (isfield(res, 'tcspcdata2')))
%         tcspcdata = sum(res.tcspcdata2(:,:,UsedBunchs),3);
%     else
%         tcspcdata = res.tcspcdata;
%     end
    if ( exist(tcspcfile, 'file') == 2) % check if tcspc data was saved
        load(tcspcfile, '-mat');
    else
         
        tst = 1;
        cnt = 0;
        bin = 0:NChannels-1;
        overcount=0;
        dind = [];
        rev_dind = zeros(8,1);
        while tst>0
             if strcmp(name(end-2:end),'ht3')
                 [tmpy, tmpx, flv, special, tst, head] = HT3_Read(name, [cnt+1, photons]);
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
                       if (dd > 0) && tmpx(j) < NChannels
                          tcspcdata(tmpx(j)+1, dd) = tcspcdata(tmpx(j)+1, dd) + 1;
                       else overcount=overcount+1;
                       end
                   end
              end
        end
 end
  
    head=HT3_Read(name);
    save(tcspcfile, 'bin', 'tcspcdata');
    close(h)
end
