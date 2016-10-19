dirname = 'm:\Narain\121113\';
names = dir([dirname '*uL*.mat']);

% read in the reference curve with no flow:
load('m:\Narain\121113\0,01nMAtto655_10min.mat');
auto0 = squeeze(sum(res.auto(:,1,2,:)+res.auto(:,2,1,:),4)); 
autotime = res.autotime;

% fit the reference curve with no-flow single-focus model:
prig = [1e-2 0.2]; % initial guess
close; 
prig = Simplex('Rigler', prig, [], [], [], [], autotime, auto0, [], 1);
% prig(1) = lateral diffusion time
% prig(2) = ratio between axial and lateral diffusion time

% read in all other curves:
j=1;
load([dirname names(j).name]);
auto = squeeze(sum(res.auto(:,1,2,:)+res.auto(:,2,1,:),4));
pos = strfind(names(j).name,'uL');
flow(j) = str2num(names(j).name(1:pos(1)-1));
for j=2:length(names)
    load([dirname names(j).name]);
    auto(:,j) = squeeze(sum(res.auto(1:size(auto,1),1,2,:)+res.auto(1:size(auto,1),2,1,:),4)); % autocorrelation curve
    pos = strfind(names(j).name,'uL');
    flow(j) = str2num(names(j).name(1:pos(1)-1)); % flow velocity
end

% fit last flow curve by using diffusion values 'prig' from data 'auto0'
close; pflow = Simplex('RiglerFlow',[prig' 100],[prig' 0],[prig' inf],[],[],autotime,auto(:,end),[],1);