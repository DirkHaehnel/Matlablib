%close all

% job = createJob(jm);
% 
% createTask(job, @sum, 1, {[1 1]});
%  
% 
% createTask(job, @sum, 1, {[2 2]});
%  
% 
% createTask(job, @sum, 1, {[3 3]});
%  
% 
% submit(job)
%  
% 
% waitForState(job, 'finished', 60)
%  
% 
% results = getAllOutputArguments(job)


jm = findResource('scheduler','configuration','hal9001')
pj = createParallelJob(jm);
createTask(pj, @labindex, 1, {});
set(pj, 'MaximumNumberOfWorkers', 3);
set(pj, 'MinimumNumberOfWorkers', 3);
submit(pj)
waitForState(pj, 'finished', 60)
results = getAllOutputArguments(pj)