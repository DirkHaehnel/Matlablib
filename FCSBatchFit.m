function [dc w0 a0 triplet c og ow oa] = FCSBatchFit(dirname, fnames, conc)

global pd
    
p = [400 150];
for j=1:length(fnames)
    eval(['load ' '''' dirname fnames(j).name '''' ' res']);
    [y, t] = FCSCrossRead(res,[1e-5 1]);
    bootstrap = [max([10 size(y,3)]) floor(size(y,3)/2)];
    pd = 1e-8./1e-6; expflag=0; bounds = [];
    data.t = t;
    data.y = y(:,1:6,:);
    [~, wtst, atst] = FCSFit(data,[400 150],expflag);
    [dc{j} w0{j} a0{j} triplet{j} c{j}] = FCSFit(data,[wtst atst],expflag,bootstrap,[],bounds);
    subplot(4,1,1:3)
    title(['Ca^{2+} conc. = ',num2str(conc(j)), ' M'])
    eval(['print -dpng -r300 ' '''' dirname fnames(j).name(1:end-4) '''' '_FCSFitBootstrap']);
    [y, t] = FCSCrossRead(res,[1e-6 1]);
    data.t = t;
    data.y = y(:,7:12,:);
    pd = 0.1; expflag=0; bounds = [];
    [og{j} ow{j} oa{j}] = FCSFit(data,[400 150],expflag,bootstrap,[],bounds);
    subplot(4,1,1:3)
    title(['Ca^{2+} conc. = ',num2str(conc(j)), ' M'])
    eval(['print -dpng -r300 ' '''' dirname fnames(j).name(1:end-4) '''' '_OregonGreenFCSFitBootstrap']);
    eval(['save ' '''' dirname '''' 'FCSFitResultsBootstrap fnames conc dc w0 a0 triplet c og ow oa']);
end

for j=1:length(dc) 
    cam(j) = mean(Diff2Rad(dc{j},20));
    cas(j) = std(Diff2Rad(dc{j},20));
    ogm(j) = mean(Diff2Rad(og{j},20));
    ogs(j) = std(Diff2Rad(og{j},20));
end
close all
% semilogx(conc,ogm,'o'); hold on; for j=1:length(ogm) line([1 1]*conc(j),[ogm(j)-ogs(j) ogm(j)+ogs(j)]); end; hold off
% figure
ind = cas<0.5;
semilogx(conc(ind),cam(ind),'o'); hold on; for j=1:length(cam) if ind(j) line([1 1]*conc(j),[cam(j)-cas(j) cam(j)+cas(j)]); end; end; hold off
xlabel({'','conc. (M)'},'interpreter','latex','fontsize',14); ylabel({'radius (\AA)',''},'interpreter','latex','fontsize',14)

