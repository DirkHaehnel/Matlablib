function [res, head] = SingleFocusCW2FCS_Cluster(name, maxtime, photons)

if nargin<5 || isempty(photons)
    photons = 1e6; % number of photons processed one at a time
end

head = ht3v2read_head(name);

%jm = findResource('scheduler','type','jobmanager','Name','hal9001','LookupURL','134.76.92.49');
jm=parcluster('SKYNET');
if (isempty(jm))
    fprintf('\n\t\tJobmanager not found\n');
    return
end

if ((isempty(head)) && not(isempty(strfind(name, '.ht3'))))
    % no header => old ht3 File
    fprintf(1,'\n\n      Warning: No File Header. Assuming old ht3 File\n');
    head.Ident = 'HydraHarp';
    head.Resolution = 2e-3;     % in nanoseconds
    head.SyncRate = 80e6;       % 80 MHz sync Rate
    head.DataStart = 0;
    tmp = dir(name);
    head.Records = tmp.bytes/4;
end;
    
Timeunit = head.Resolution;                 % timeunit in ns
Timeunit2 = 1/(head.SyncRate);              % timeunit in s
microresolution = head.Resolution;          % resolution in ns
syncperiod = (1 / head.SyncRate) * 1e9;     % syncperiod in ns

if nargin<2 || isempty(maxtime)
    maxtime = 1e3;
    Nsub = 10;
    Ncasc = ceil(log2(1e3/Timeunit/Nsub));    % 1 µsec max correlation time
    Ncasc2 = ceil(log2(3/Timeunit2/Nsub));    % 3 sec max correlation time
else
    if length(maxtime)==2
        Nsub = maxtime(2);
    else
        Nsub = 10;
    end
    Ncasc = ceil(log2(maxtime(1)/Timeunit/Nsub));
    Ncasc2 = ceil(log2(3/Timeunit2/Nsub));    % 3 sec max correlation time
end
autotime = cumsum(reshape(repmat(2.^(0:Ncasc-1),Nsub,1),Ncasc*Nsub,1));
autotime = autotime*Timeunit;      % in ns

autotime2 = cumsum(reshape(repmat(2.^(0:Ncasc2-1),Nsub,1),Ncasc2*Nsub,1));
autotime2 = autotime2*Timeunit2;   % in s

autowidth = sum(horzcat(vertcat(0, diff(autotime)/2), vertcat(diff(autotime)/2, 0)),2); % in ns
autowidth = autowidth * 1e-9;   % in s

%autotime = 0:160:ceil(automaxtime*1e3);  % in ps
%autotime = autotime*1e-3;                % in ns

NCounts = head.Records;

[sync, micro, chan, special] = ht3v2read_old(name, [1, photons]);
tchan = chan(special == 0);
dind = unique(tchan(1:1e4));
ndet = length(dind);
rev_dind = zeros(8,1);
for j=1:length(dind)
    rev_dind(dind(j)+1) = j;
end

tcspc = zeros(2^15, ndet);
jj = 1;
cnt = 0;
tst = 1;

% Randompart for the filename
random_name = round(rand * 10000000);

% number of bunches
anzbunch = ceil( NCounts / photons);

time = zeros(anzbunch, 1);
rate = zeros(anzbunch, ndet);
auto = 0*autotime;
auto = repmat(auto,[1 ndet ndet anzbunch]);
auto2 = zeros(length(autotime2), ndet, ndet, anzbunch);

jobs = cell(anzbunch, 1);
filelist = cell(anzbunch, 1);

h = waitbar(0,'Starting Jobs...');
while (tst>0)
    [sync, micro, chan, special, tst] = ht3v2read_old(name, [cnt+1, photons]);
    if (tst>0) && (~isempty(sync))
        cnt = cnt + tst;
        
        % save for cluster
        temp_name = [name(1:end-4) sprintf('_%i-%i.mat', random_name, jj)];
        s2 = regexp(temp_name, '\', 'split');
        file_name = s2(end);
        save(temp_name, 'sync', 'micro', 'chan', 'special', 'syncperiod', 'microresolution', 'dind', 'rev_dind', 'maxtime', 'ndet', 'autotime', 'autowidth', 'autotime2', 'Timeunit2', 'Ncasc2', 'Nsub');
        filelist{jj} = temp_name;
        
        filedeps = cell(1,4);
        filedeps{1} = 'tttr2xfcs.m';
        filedeps{2} = 'SingleFocusCW2FCS_Cluster_Calc.m';
        filedeps{3} = 'mHist.m';       
        filedeps{4} = temp_name;
        
        jobs{jj} = createJob(jm);
        % We missuse the JobData Property to store the location for
        % the result
        set(jobs{jj}, 'JobData', jj);
        set(jobs{jj}, 'FileDependencies', filedeps);
        
        createTask(jobs{jj}, @SingleFocusCW2FCS_Cluster_Calc, 5, { file_name });

        if (length(findJob(jm, 'State', 'queued')) < 20)
            submit(jobs{jj});
        end
        waitbar(jj/anzbunch,h)
        jj = jj+1;
    end
end
submittedJobs = jj - 1;
close(h)

finished_jobs = 0;
h = waitbar(0,'Finished Jobs...');
while (finished_jobs < submittedJobs)
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
    
    waitbar(finished_jobs/submittedJobs,h)
    fprintf('%i Jobs of %i finished\n', finished_jobs, length(jobs));
    if (finished_jobs < submittedJobs)
       pause(10) 
    end    
end
close(h)

% Fetch the data - when all jobs have finished
for jj=1:submittedJobs
    %jobdata = get(jobs{jj}, 'JobData');
    outargs = getAllOutputArguments(jobs{jj});

    auto(:,:,:,jj) = outargs{1};
    auto2(:,:,:,jj) = outargs{2};
    tcspc = tcspc + outargs{3};
    rate(jj,:) = outargs{4};
    time(jj) = outargs{5};
    
    destroy(jobs{jj});
end

% Finally remove the intermediatly created .mat files
for i=1:length(filelist)
    delete(filelist{i});
end

clf
blaa.auto = auto;
blaa.autotime = autotime*1e-9;
blaa2.auto = auto2;
blaa2.autotime = autotime2;
[y, t] = FCSCrossReadCW(blaa);
[y2, t2] = FCSCrossReadCW(blaa2);
semilogx(t, mean(y,3), t2(3:end), mean(y2(3:end,:,:),3));
%semilogx(autotime*1e-9, reshape(sum(auto,4),[size(auto,1) size(auto,2)*size(auto,3)])/sum(time));

tau = (0:2^15-1)*microresolution;
tau = tau';
tau(any(tcspc==0,2),:)=[];
tcspc(any(tcspc==0,2),:)=[];

res.tau = tau;
res.tcspc = tcspc;

res.autotime1 = autotime;
res.autotime2 = autotime2;
res.auto1 = auto;
res.auto2 = auto2;

[res.autotime ind] = sort(vertcat(autotime*1e-9, autotime2(3:end)));   % in s
tmp = cat(1, auto(:,:,:,:), auto2(3:end,:,:,:));
res.auto = tmp(ind,:,:,:);

res.autowidth1 = autowidth;
res.autowidth = sum(horzcat(vertcat(0, diff(res.autotime)/2), vertcat(diff(res.autotime)/2, 0)),2); % in s

res.rate = rate;
res.time = time;

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



