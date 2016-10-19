function[autotime auto automean autofinal master] = ht3FLCS_Pro(irf, name)
% This function genetaes Filters for the fitted lifetime curves to be given
% as input for Fluorescence Lifetime Correlation Spectroscopy
% name is the file input for the TCSPC decay curve
% param is a matrix containing fitted lifetime curves 
% program works for 2 detectors
% p gives the max TCSPC time used to record the data and dt is the
% resolution.

% (c) Narain Karedla 2013

photons=1e6; % 1 million photons at a time
%timeframe=60;
%% Obtaining IRF data
if strcmp(irf(end-2:end), 'ht3')
    headirf=HT3_Read(irf);
    Timeunitirf=1/(headirf.SyncRate);
    Resolutionirf=headirf.Resolution*1e-9;
    NChannelsirf=ceil(Timeunitirf/Resolutionirf);
    NCountsirf=headirf.Records;
else
    disp('Not a valid file!')
    return
end
       % h = waitbar(0,'Please wait...');
        tcspcirf = [];
        tcspcfile = [irf(1:end-7) '.ht3tcspc'];
        matfile = [irf(1:end-4) '.mat'];

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
        tcspcdata=[];
        tst = 1;
        cnt = 0;
        tau = 0:NChannelsirf-1;
        overcount=0;
        dind = [];
        rev_dind = zeros(8,1);
        while tst>0
             if strcmp(irf(end-2:end),'ht3')
                 [tmpy, tmpx, flv, special, tst, head] = HT3_Read(irf, [cnt+1, photons]);
             end
             if isempty(dind)
                 %dind = union(flv,flv);
                 dind = unique(flv);
                 dnum = length(dind);
                 
                     tcspcdata = zeros(NChannelsirf,dnum);
                 
                    for j=1:length(dind)
                         rev_dind(dind(j)+1) = j;
                    end
             end
                 cnt = cnt+tst;
                % waitbar(cnt/NCountsirf,h)
              if (tst>0) && (~isempty(tmpy))
                   for j=1:length(flv)
                        dd = rev_dind(flv(j)+1);
                       if (dd > 0) && tmpx(j) < NChannelsirf
                          tcspcdata(tmpx(j)+1, dd) = tcspcdata(tmpx(j)+1, dd) + 1;
                       else overcount=overcount+1;
                       end
                   end
              end
        
       end
   % close(h)
    save(tcspcfile, 'tau', 'tcspcdata');
    end
    tcspcirf=tcspcdata;
    tcspcdata=[];

%%    Obtaining TCSPC data
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
    %h = waitbar(0,'Please wait...');
    tau = [];
    tcspcdata = [];
    tcspcfile = [name(1:end-7) '.ht3tcspc'];
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
         
        tst = 1;
        cnt = 0;
        tau = 0:NChannels-1;
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
                 %waitbar(cnt/NCounts,h)
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
    %close(h)
    save(tcspcfile, 'tau', 'tcspcdata');
    end

   
%% Lifetime fitting
    if Resolutionirf ~=Resolution && sigdisp(Timeunitmain,4) ~= sigdisp(Timeunit,4)
        disp ('Resolution mismatch!! Use correct irf file!')
    return
    end
    
    if size(tcspcdata,2)~=size(tcspcirf,2)
        disp('Detectors mismatch')
        return
    end
    
     % deleting extra channels at the end
            ind=1:size(tcspcirf,1);
            tmp=ind>NChannelsirf;
            tcspcirf(tmp==1,:)=[];
            
            ind=1:size(tcspcdata,1);
            tmp=ind>NChannels;
            tcspcdata(tmp==1,:)=[];
            tau=((0:NChannels-1).*Resolution)';
            if size(tcspcdata,1)>size(tcspcirf,1)
                tcspc(size(tcspcirf,1)+1:end,:)=[];
            else
                tcspcirf(size(tcspcdata,1)+1:end,:)=[]; 
            end 
            
        % deriving lifetime values
            p=ceil(Timeunit*1e9);
            dt=Resolution*1e9;
            k=size(tcspcirf,2);
            l=size(tcspcirf,1);
            
            for i=1:k;
                    [~, offset(:,i), ~, life(:,1), ~, dlife(:,1), irs(:,1), zz(:,:,1), ~, ~] = Fluofit(abs(tcspcirf(1:end,i)-mean(tcspcirf(end-1000:end-10,i))), tcspcdata(1:end,i), p,dt);     
                     mattau{i}=life;
                     matdtau{i}=dlife;
                     matzz{i}=zz;
                     matirs{i}=irs;
                     
        
                     life=[]; zz=[]; dlife=[]; irs=[]; 
            end
%            for i=1:size(matzz,2)
%                for j=1:size(matzz{i},2)
%                if mean(matzz{i}(:,j))<0;
%                    matzz{i}(:,j)=[];
%                end
%            end
           
              c=0; d=0; 
              param=[];
             for i=1:size(matzz,2)
                 for j=1:size(matzz{i},2)
                    c=c+1;
                    param(:,c)=matzz{i}(:,j);
                    %irs(:,c)=matirs{i}(:,j);
                 end
                 for e=1:size(mattau{i})
                     d=d+1;
                     life(d,:)=[mattau{i}(e) i];
                     dlife(d,:)=[matdtau{i}(e) i];
                     
                 end 
             end
             %o=0;
             for i=1:size(param,2)
                 if mean(param(:,i))<1
                     param(:,i)=1;
                     %o=o+1
                 end
             end
             
    %% Performing FLCS
    for m=1:size(param,2)
        param(:,m) = param(:,m)/sum(param(:,m));
    end % normalization
    c=0;
    master=[];
    det=size(tcspcdata,2);
    k=size(param,2)/det;    

    for i=1:det;
        c=c+1;
        l=1+(c-1)*k;
        master = [master, ((param(:,l:c*k)'*diag(1./tcspcdata(:,c))*param(:,l:c*k))\(diag(1./tcspcdata(:,c))*param(:,l:c*k))')'];
    end
    time = 0;
    cnt = 0;
    tst = 1; 
    jj=1;
    flv=[]; tmpx=[]; tmpy=[]; 
    Nsub  = 10;
    Ncasc    = ceil(log2(3/Timeunit/Nsub+1)); % 3 sec max correlation time
    autotime = cumsum(reshape(repmat(2.^(0:Ncasc-1),Nsub,1),Ncasc*Nsub,1));
    autotime = autotime*Timeunit;
    auto=zeros(length(autotime), size(master,2), size(master,2));
    
    close all;
        
    while tst>0
        
                [tmpy, tmpx, flv, markers, tst, loc, head1] = HT3_Read(name,[cnt+1 photons]); 
                if (tst>0)
                ind= flv>=det;
                tmpy(ind) = [];
                tmpx(ind) = [];
                flv(ind)  = [];
                
                cnt=cnt+tst;
                if tst>0 
                    dt=(tmpy(end)-tmpy(1)*Timeunit)+(tmpx(end)-tmpx(1)*Resolution); 
                    time(jj)= dt*Timeunit;
                    
                    tmpx  = tmpx*Resolution;    % in seconds
                    
                    tmpy  = (tmpy*Timeunit+tmpx);  %in seconds
                    
                    tmp   = tmpx>=tau(1) & tmpx<=tau(end);
                    tmpy  = tmpy(tmp);
                    tmpx  = tmpx(tmp);
                    flv   = flv(tmp);
                    
                       ind0=(flv==0);
                       ind1=(flv==1);
                    
                    if numel(tmpy>0)
                       tmpy=(tmpy-tmpy(1))/Timeunit;  % back into sync units
                       tmpx=round((tmpx-tau(1))/Resolution)+1; % TCSPC Channel
                       
                       
                       p=size(master,2)/det;
                       q=size(flv,1);
                       flvp = [];
                       flvp(1:q,1:2*p)=0;
                       j=[];
                       for j=1:p;
                           
                           flvp(:,j)=master(tmpx,j).*ind0;
                           flvp(:,j+p)=master(tmpx,j+p).*ind1;
                           
                       end
                       
                    end
                end
                auto(:,:,:,jj)=tttr2xfcs(tmpy,flvp,Ncasc,Nsub);
                subplot('position', [0.925 0.2 0.025 0.6]);
                bar(0,cnt/NCounts);
                axis([-0.4 0.4 0 1]);
                set(gca,'xtick',[],'ytick',[]);
                subplot('position', [0.1 0.1 0.7 0.8]);
                semilogx(autotime, reshape(sum(auto,4),[size(auto,1) size(auto,2)*size(auto,3)])/sum(time));
                drawnow;
                end
                jj=jj+1;
                   
    end
        for r=1:length(time)
                        auto(:,:,:,r) = auto(:,:,:,r)/time(r);
        end
clf
semilogx(autotime, reshape(sum(auto,4),[size(auto,1) size(auto,2)*size(auto,3)])/sum(time));      
automean=mean(auto,4);

  for s=1:size(master,2)/det;
      autofinal(:,s)=(automean(:,s,s+p)+automean(:,s+p,s))/2;
      autofinal(:,s)=autofinal(:,s)-min(autofinal(end-20:end,s));
      autofinal(:,s)=autofinal(:,s)/max(autofinal(:,s))+s;
  end
    
  subplot(2,1,1)
  plot(tau, master)
  subplot(2,1,2)
  semilogx(autotime,autofinal)
end
