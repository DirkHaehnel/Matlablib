function [res, head] = MultiFocus2FCS_Cluster(name, cnum, maxtime, para, cutoff, photons)

if nargin<6 || isempty(photons)
    photons = 1e6; % number of photons processed one at a time
end

close all

%jm = findResource('scheduler','type','jobmanager','Name','hal9001','LookupURL','134.76.92.49');
jm=parcluster('SKYNET');
if (isempty(jm))
    fprintf('\n\t\tJobmanager not found\n');
    return
end

if strcmp(name(end-2:end),'pt3')
    head = pt3v2Read(name);
    Timeunit = head.Sync*1e-9;
    Resolution = head.Resolution;
    NChannels  = 2^12;
    NCounts   = head.Records;
elseif strcmp(name(end-2:end),'t3r')
    head = tttrRead(name);
    Timeunit = head.GlobClock*1e-9;
    Resolution = head.Resolution;
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
        Timeunit = 1/(head.SyncRate);
        Resolution = head.Resolution*1e-9;
        NChannels = 2^15;
        tmp = dir(name);
        NCounts = tmp.bytes/4;
    end
else
    disp('Not a valid file!');
    return
end

if nargin<3 || isempty(maxtime)
    maxtime = 3;
    Nsub = 10;
    Ncasc = ceil(log2(maxtime/Timeunit/Nsub)); % 3 sec max correlation time
else
    if length(maxtime)==2
        Nsub = maxtime(2);
    else
        Nsub = 10;
    end
    Ncasc = ceil(log2(maxtime(1)/Timeunit/Nsub));
end
%if 3e3*maxtime(1)*NCounts/head.AcquisitionTime<photons
%    photons = round(3e3*maxtime(1)*NCounts/head.AcquisitionTime);
%end
autotime = cumsum(reshape(repmat(2.^(0:Ncasc-1),Nsub,1),Ncasc*Nsub,1));
autotime = autotime*Timeunit;

if strcmp(name(end-2:end),'pt3')
    [tmpy, tmpx, flv, bla, tst] = pt3v2Read(name, [1, 1e5]);
elseif strcmp(name(end-2:end),'t3r')
    [tmpy, tmpx, flv, bla, tst] = tttrRead(name, [1, 1e5]);
elseif strcmp(name(end-2:end),'ht3')
    if isempty(head)
        [tmpy, tmpx, flv, special, tst, overcount, head] = ht3read(name, [1, 1e5]);
    else
        [tmpy, tmpx, flv, special, tst, head] = ht3v2read(name, [1, 1e5]);
    end
end

dind = unique(flv);
dnum = length(dind);
for j=1:length(dind)
    rev_dind(dind(j)+1) = j;
end
        
% Read TCSPC Data
[bin, tcspcdata] = ReadTCSPC(name);
tcspcdata2 = zeros(NChannels,dnum,ceil(head.Records/photons));

% Detect Laser Puls Time Gates
[t1, len] = AutodetectTimeGates(tcspcdata, cnum);

for j=1:cnum
    for k=1:dnum
        tau(1:len,j,k) = bin(t1(j):t1(j)+len-1);
        tcspc(1:len,j,k) = tcspcdata(t1(j):t1(j)+len-1,k);
    end
end
if nargin>4 && ~isempty(cutoff)
    ind = tau(:,1,1) - tau(1,1,1) > cutoff/Resolution;
    tau = tau(ind,:,:);
    tcspc = tcspc(ind,:,:);
end
tcspc(any(any(tau==0,2),3),:,:)=[];
tau(any(any(tau==0,2),3),:,:)=[];
semilogy(Resolution*tau(:,:,1),tcspc(:,:,1)); drawnow
res.tau = Resolution*tau;
res.tcspc = tcspc;
res.bin = bin;
res.tcspcdata = tcspcdata;

rate = [];
time = 0;
cnt = 0;
tst = 1;
jj = 1;

% Randompart for the filename
random_name = round(rand * 10000000);

% number of bunches
anzbunch = ceil( NCounts / photons);

time = zeros(1, anzbunch);
rate = zeros(anzbunch, dnum*cnum);
auto = zeros(length(autotime), cnum*dnum, cnum*dnum, anzbunch);
jobs = cell(anzbunch, 1);
filelist = cell(anzbunch, 1);

h = waitbar(0,'Starting Jobs...');
while tst>0
    if strcmp(name(end-2:end),'pt3')
        [tmpy, tmpx, flv, bla, tst] = pt3v2Read(name, [cnt+1, photons]);
    elseif strcmp(name(end-2:end),'t3r')
        [tmpy, tmpx, flv, bla, tst] = tttrRead(name, [cnt+1, photons]);
    elseif strcmp(name(end-2:end),'ht3')
        if isempty(head)
            [tmpy, tmpx, flv, bla, tst] = ht3read(name, [cnt+1, photons]);
        else
            [tmpy, tmpx, flv, bla, tst] = ht3v2read(name, [cnt+1, photons]);
        end
    end
    cnt = cnt + tst;
    if (tst>0) && (~isempty(tmpy))
        dt = tmpy(end)-tmpy(1);
        time(jj) = dt*Timeunit;
        
        % save for cluster
        temp_name = [name(1:end-4) sprintf('_%i-%i.mat', random_name, jj)];
        s2 = regexp(temp_name, '\', 'split');
        file_name = s2(end);
        save(temp_name, 'tmpy', 'tmpx', 'flv', 'tau', 'dind', 'rev_dind', 'dnum', 'cnum', 'Ncasc', 'Nsub', 'NChannels', 'dt', 'Timeunit');
        filelist{jj} = temp_name;
        
        filedeps = cell(1,3);
        filedeps{1} = 'tttr2xfcs.m';
        filedeps{2} = 'MultiFocus2FCS_Cluster_Calc.m';
        filedeps{3} = temp_name;
        
        jobs{jj} = createJob(jm);
        % We missuse the JobData Property to store the location for
        % the result
        set(jobs{jj}, 'JobData', jj);
        set(jobs{jj}, 'FileDependencies', filedeps);
        
        createTask(jobs{jj}, @MultiFocus2FCS_Cluster_Calc, 3, { file_name });

        if (length(findJob(jm, 'State', 'queued')) < 20)
            submit(jobs{jj});
        end
        waitbar(jj/anzbunch,h)
    end
    jj = jj+1;
end
close(h)

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
            [pending running finished] = findTask(jobs{jj});
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
    
    waitbar(finished_jobs/anzbunch,h)
    fprintf('%i Jobs of %i finished\n', finished_jobs, length(jobs));
    if (finished_jobs < length(jobs))
       pause(10) 
    end    
end
close(h)

% Fetch the data - when all jobs have finished
for jj=1:length(jobs)
    %jobdata = get(jobs{jj}, 'JobData');
    outargs = getAllOutputArguments(jobs{jj});

    tcspcdata2(:,:,jj) = outargs{3};
    rate(jj,:) = outargs{2};
    auto(:,:,:,jj) = outargs{1};
    
    destroy(jobs{jj});
end

% Finally remove the intermediatly created .mat files
for i=1:length(filelist)
    delete(filelist{i});
end

for j=1:length(time)
    auto(:,:,:,j) = auto(:,:,:,j)/time(j)/Timeunit;
end

clf
semilogx(autotime, reshape(sum(auto,4),[size(auto,1) size(auto,2)*size(auto,3)])/sum(time));

ind = sum(sum(tcspcdata2,3),2)==0;
tcspcdata2(ind,:,:) = [];
res.tcspcdata2 = tcspcdata2;

res.autotime = autotime;
res.auto = auto;
res.rate = rate;
res.time = time;
head.NCounts = cnt;
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
