function [dc w0 a0 triplet c velo err z] = FCSFit_Cluster(name,p,expflag,bs,para,bounds,pmin,pmax,bld, timebounds)

% Function FCSFit for fitting 2fFCS data
% INPUT PARAMETERS
% name       : Is the name of a file on the hard disk containing a structure (res) as
%              generated by programs TwoFocus2FCS.m or Multifocus2FCS.m. 
% p          : Initial guess value vector. p(1) for beam waist radius. p(2) for
%              confocal pinhole radius. p(3:5) are the flow velocity parameters.
%              Default is p=[400 150]
% expflag    : For fitting exponential decay terms. The routine assumes that
%              pd(end-expflag+1:end) are the initial guess values of exponential 
%              decay times in microseconds whereas the other pd(1:end-expflag) 
%              values are guess values for the inverse of diffusion coefficients.
% bs         : If bs=1, the program does bootstrapping in the cluster. If
%              bs=0, the sum of the bunches are taken together for fitting.
% para       : The Optical setup vector: para(1) is radius of pinhole divided by
%              the magnification. para(2) is excitation wavelength divided by the 
%              refractive index of the immersion medium. para(3) is center
%              emission divided by the refractive index of the immersion
%              medium. para(4) is the interfocal distance. All in nanometers.
% bounds     : bounds on the fitting variable pd. bounds must be 2Xlength(pd).
% pmin       : lower bounds of p during fitting.
% pmax       : upper bounds of p during fitting.
% timebounds : time range for the data to be fitted. By default, all data 
%              from 0 to 230 is fitted.  
% 
% OUTPUT PARAMETERS
% dc      : diffusion coefficients in sqcm/sec.
% w0      : beam waist radius.
% a0      : confocal pinhole aperture.
% triplet : exponential decay times of triplet states.
% c       : amplitudes of different components contributing to correlation.
% velo    : velocity vector
% err     : fit error
% z       : fitted curve

% (c) Narain Karedla 2013

close all
global pd
 
%pd=0.1; 


random_name = round(rand * 10000000);

if nargin<2 || isempty(p)
    p = [400 150];
end

if nargin>4 && ~isempty(para)
    if ischar(para)
        if strcmp(para,'b')
            av = 75e3/60;
            lam = [470 530]/1.33;
            delta = 393;
        elseif strcmp(para,'r')
            av = 75e3/60;
            lam = [640 670]/1.33;
            delta = 414;
        else
            disp('wrong parameter')
            return;
        end
    else
        av = para(1);
        lam = para(2:3);
        delta = para(4);
    end
else
    av = 75e3/60;
    lam = [640 670]/1.33;
    delta = 414;
end

if length(p)==3 && p(3)==1i
    p = p(1:2);
    flag2d = 1;
else
    flag2d = 0;
    if length(p)>2
        p(3:end) = pd(1)*p(3:end)/1e3;
    end
end

% if isstruct(name)
%     if isfield(name,'auto')
%         [y,t] = FCSCrossRead(name);
%     else
%         y = name.y;
%         t = name.t;
%     end
% else
%     [y,t] = FCSCrossRead(name);
% end

load(name);
[y,t]= FCSCrossRead(res);
if ~isempty(timebounds)
    y=y(timebounds(1):timebounds(2),:,:);
    t=t(timebounds(1):timebounds(2));
end

   
if size(y,3)>1
    y = y(:,:,~isnan(sum(sum(y,1),2)));
end

if isempty(pd)
    disp('Please define global variable pd');
    return
end
if nargin<3 || isempty(expflag)
    expflag = [length(pd)-1 0];
elseif length(expflag)==1
    expflag = [expflag 0];
end
if nargin<6 || isempty(bounds)
    bounds = [];
end
if nargin<7 || isempty(pmin)
    pmin = 0*p;
    if length(p)>2
        pmin(3:end) = -inf;
    end
end
if nargin<8 || isempty(pmax)
    pmax = inf*ones(size(p));
end
if nargin<9 
    bld = 1;
end

%% bootstrapping...
if bs==1
 %   jm = findResource('scheduler','type','jobmanager','Name','hal9001','LookupURL','134.76.92.49');
jm=parcluster('SKYNET');
    if (isempty(jm))
        fprintf('\n\t\tJobmanager not found\n');
        return
    end
bootstrap=[ceil(size(y,3)/2) round(size(y,3)/2)];
    if size(y,3)>1
        weight = 1./std(y,0,3);
    else
        weight = [];
    end
    p0 = Simplex('GaussFcs',p,pmin,pmax,[],[],av,lam,delta,1e6*t,mean(y(:,1:2,:),3),mean(y(:,3:end,:),3),expflag,bounds,flag2d,weight,bld);
    pd0 = pd;
    
    % bootstrap begins
    for j=1:bootstrap(1)
        pd=pd0;
        % real bootstrap choosing bootstrap(1) times bootstrap(2) random data
        ind = randi(size(y,3),[1 bootstrap(2)]);
        
        databs.y=mean(y(:,:,ind),3);
        databs.t=t;
        
        temp_name = [name(1:end-4) sprintf('_%i-%i.mat', random_name, j)];
        s2 = regexp(temp_name, '\', 'split'); %cell array of D:\, folder1, folder2, ... file name
        file_name = s2(end);
        save(temp_name,'databs','av','lam','bld','bounds','delta','expflag','flag2d','lam','p0', 'pd0','pmin','pmax','weight');
        filelist{j} = temp_name;
        
        filedeps = cell(1,7);
        filedeps{1} = 'Simplex.m';
        filedeps{2} = 'GaussFcs.m';
        filedeps{3} = 'FCSFit_Cluster_Calc.m';
        filedeps{4} = 'LaserBeamFit.m';
        filedeps{5} = 'D1CEF.m';
        filedeps{6} = 'AffineFit.m';
        filedeps{7} = temp_name;
        
        jobs{j} = createJob(jm);
        set(jobs{j}, 'JobData', j);
        set(jobs{j}, 'FileDependencies', filedeps);
        createTask(jobs{j}, @FCSFit_Cluster_Calc, 9, { file_name });
        
        
        if (length(findJob(jm, 'State', 'queued')) < 20)
            submit(jobs{j});
        end
                
    end
    
    finished_jobs = 0;
    h = waitbar(0,'Finished Jobs...');
    while (finished_jobs < length(jobs))  
        finished_jobs = 0;
       for jj=1:length(jobs)
        if strcmp(get(jobs{jj}, 'State'), 'failed')
            % Seems to be necessary as the RZ Aachen has the brain-dead
            % config of running more workers than licences available
            jobs{jj} = resubmitjob(jm, jobs{jj});
            continue;
        elseif strcmp(get(jobs{jj}, 'State'), 'finished')
            [~,~,finished] = findTask(jobs{jj});
            without_error = 1;
            for j=1:size(finished, 2)
                if ~ isempty(get(finished(j), 'ErrorIdentifier'))
                    fprintf('Will resubmit job %i because of an Error in Task %i:\n%s\n', jj, j, get(finished(j), 'ErrorMessage'));
                    without_error = 0;
                end
            end
            if without_error == 1 || size(finished, 2) == 0
                finished_jobs = finished_jobs + 1;
            else
                jobs{jj} = resubmitjob(jm, jobs{jj});
                continue;
            end
        elseif strcmp(get(jobs{jj}, 'State'), 'pending')
            if (length(findJob(jm, 'State', 'queued')) < 20)
                submit(jobs{jj});
            end
        end
       end
        waitbar(finished_jobs/bootstrap(1),h)
        fprintf('%i Jobs of %i finished\n', finished_jobs, length(jobs));
        if (finished_jobs < length(jobs))
         pause(60) %pause for jobs to finish
        end    
    end
    close(h);
    
    p=[];
    for j=1:length(jobs)
   
        outargs = getAllOutputArguments(jobs{j});
        p(:,j)=outargs{1};
        err(:,j)=outargs{2};
        c(:,:,j)=outargs{3};
        if expflag~=0;
            triplet(:,j)=outargs{4};
        else
            triplet=[];
        end
        velo(:,j)=outargs{5};
        w0(:,j)=outargs{6};
        a0(:,j)=outargs{7};
        dc(:,j)=outargs{8};
        z(:,:,j)=outargs{9};
    
        destroy(jobs{j});
    end
    
    for i=1:length(filelist)
     delete(filelist{i});
    end
    
     z_av=mean(z,3); % Average z fit of all the jobs
     y_av=mean(y,3); % Average of all the correlation data
     y_av(:,3)=y_av(:,5)+ y_av(:,3);
     y_av(:,4)=y_av(:,6)+ y_av(:,4);
     y_av(:,5:6)=[];
%     for i=1:2;
%         p0(i)=mean(p(i,:));
%     end
   
     subplot(4,1,1:3)
     semilogx(t,y_av,'o', t,z_av)
    %[tmp, tmp, z] = GaussFcs(p0,av,lam,delta,1e6*t,mean(y(:,1:2,:),3),mean(y(:,3:end,:),3),expflag,bounds,flag2d,mean(weight,3));
    
     ax=axis;axis([min(t) max(t) ax(3:4)]);
      ylabel('correlation / (\ith\rm\nu)^2/s^2'); 
     for j=1:length(pd)-sum(expflag)
        tmp = -floor(log10(mean(dc(j,:))));
        s{j} = ['\itD\rm_' mint2str(j) ' = (' mnum2str(mean(dc(j,:))*10^tmp,1,2)  '\pm' mnum2str(std(dc(j,:))*10^tmp,1,2) ')\cdot10^{-' int2str(tmp) '} cm^2/s'];
    end
    if ~isempty(velo)
        if ~(mean(velo(1,:))==0)
            s{end+1} = ['\itv_x\rm = (' mnum2str(mean(velo(1,:)),1,2) '\pm' mnum2str(std(velo(1,:)),1,2) ') \mum/s'];
        end
        if ~(mean(velo(2,:))==0)
            s{end+1} = ['\itv_y\rm = (' mnum2str(mean(velo(2,:)),1,2) '\pm' mnum2str(std(velo(2,:)),1,2) ') \mum/s'];
        end
        if ~(mean(velo(3,:))==0)
            s{end+1} = ['\itv_z\rm = (' mnum2str(mean(velo(3,:)),1,2) '\pm' mnum2str(std(velo(3,:)),1,2) ') \mum/s'];
        end
    end
    if ~isempty(triplet)
        s{end+1} = '';
        if size(triplet,1)==1
            s{end+1} = ['\tau = ' int2str(round(mean(triplet))) ' \mus'];
        else
            tmp = ['\tau = (' int2str(round(mean(triplet(1,:))))];
            for k=2:size(triplet,1)
                tmp = [tmp ',' int2str(round(mean(triplet(k,:))))];
            end
            tmp = [tmp  ') \mus'];
            s{end+1} = tmp;
        end
    end
    s{end+1} = '';
    s{end+1} = ['\itw\rm_0 = ' int2str(round(mean(w0))) ' nm'];
    if a0>0
        s{end+1} = ['\ita\rm_0 = ' int2str(round(mean(a0))) ' nm'];
    end
    text(0.6,0.95,s,'VerticalAlignment','Top','units','normal')
    subplot(4,1,4)
%     if size(y,2)==6
%         y(:,3,:)=sum(y(:,[3 5],:),2);
%         y(:,4,:)=sum(y(:,[4 6],:),2);
%         y(:,5:6,:) = [];
%     end
    semilogx(t,y_av-z_av)
    ax=axis;axis([min(t) max(t) -max(abs(ax(3:4))) max(abs(ax(3:4)))]); xlabel('time [s]'); ylabel('residuals');
    set(gcf,'Position',[20 50 900 850],'PaperPositionMode','auto')
    drawnow
    
    
elseif bs==0 
    %% Sum all fitting begins
    if size(y,3)>1
        weight = 1./std(y,0,3);
    else
        weight = [];
    end
    y = mean(y,3);
    p = Simplex('GaussFcs',p,pmin,pmax,[],[],av,lam,dist,1e6*t,y(:,1:2),y(:,3:end),expflag,bounds,flag2d,weight,bld);
    subplot(4,1,1:3)
    [err, c, z] = GaussFcs(p,av,lam,dist,1e6*t,y(:,1:2),y(:,3:end),expflag,bounds,flag2d,weight);
    if length(p)>2
        velo = zeros(3,1);
        velo(1:length(p)-2) = p(3:end)/pd(1)*1e3;
        p(3:end) = [];
    else
        velo = [];
    end
    if flag2d
        w0(1) = p(1);
        w0(2) = p(2);
        p(2) = sqrt(sum(p(1:2).^2)/2); p(1) = [];
    else
        w0 = p(1);
    end
    if length(p)>1
        a0 = p(2);
    else
        a0 = 0;
    end
    dc = 1e-14/1e-6./pd(1:end-sum(expflag));
    if expflag(2)>0
        dc(end+1:end+expflag(2)) = pd(end-sum(expflag)+1:end-expflag(1));
    end
    ax=axis;axis([1e6*min(t) 1e6*max(t) ax(3:4)]); ylabel('correlation / (\ith\rm\nu)^2/s^2'); 
    if expflag(1)>0
        triplet = pd(end-expflag+1:end);
    else
        triplet = [];
    end
    for j=1:length(pd)-sum(expflag)
        tmp = -floor(log10(dc(j)));
        s{j} = ['\itD\rm_' mint2str(j) ' = ' mnum2str(dc(j)*10^tmp,1,2) '\cdot10^{-' int2str(tmp) '} cm^2/s'];
    end
    if ~isempty(velo)
        s{end+1} = ['\bfv\rm = (' mnum2str(velo(1),1,2) ', ' mnum2str(velo(2),1,2) ', ' mnum2str(velo(3),1,2) ') \mum/s'];
    end
    if ~isempty(triplet)
        s{end+1} = '';
        if length(triplet)==1
            s{end+1} = ['\tau = ' int2str(round(triplet(1))) ' \mus'];
        else
            tmp = ['\tau = (' int2str(round(triplet(1)))];
            for k=2:length(triplet)
                tmp = [tmp ',' int2str(round(triplet(k)))];
            end
            tmp = [tmp  ') \mus'];
            s{end+1} = tmp;
        end
    end
    s{end+1} = '';
    s{end+1} = ['\itw\rm_0 = ' int2str(round(w0)) ' nm'];
    if a0>0
        s{end+1} = ['\ita\rm_0 = ' int2str(round(a0)) ' nm'];
    end
    text(0.6,0.95,s,'VerticalAlignment','Top','units','normal')    
    subplot(4,1,4)
    if size(y,2)==6
        y(:,3)=sum(y(:,[3 5]),2);
        y(:,4)=sum(y(:,[4 6]),2);
        y(:,5:6)=[];
    end
    semilogx(1e6*t,y-z)
    ax=axis;axis([1e6*min(t) 1e6*max(t) -max(abs(ax(3:4))) max(abs(ax(3:4)))]); xlabel('time [\mus]'); ylabel('residuals');
    pos = get(0,'ScreenSize');
    set(gcf,'Position',[0.1*pos(3) 0.05*pos(4) 0.65*pos(3) 0.85*pos(4)],'PaperPositionMode','auto')
    drawnow
else
    disp('The program takes two inputs for bs (1 and 0) for bootstrapping or sum fitting respectively')
    return
end
   
end




function [job] = resubmitjob(jm, old_job)
    fprintf('BEWARE -- WE RESUBMIT A JOB - THIS SHOULD NOT HAPPEN!\n');
    fprintf('YOU HAVE BEEN WARNED!\n');
    % Recreate the job state of the old_job in the new job
    % This is plain evil and against everything good!
    job = createJob(jm);
    set(job, 'JobData', get(old_job, 'JobData'));
    set(job, 'FileDependencies', get(old_job, 'FileDependencies'));
    old_tasks = get(old_job, 'Tasks');
    for i=1:size(old_tasks)
        old_task = old_tasks(i);
        createTask(job, get(old_task, 'Function'), get(old_task, 'NumberOfOutputArguments'), get(old_task, 'InputArguments'));
    end
    destroy(old_job);
    submit(job);
end
