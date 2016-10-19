% AntibunchAnalysis1

names = {'GFPl','2GFPl','rGFPtriplo','hGFPtril','htriGFPl','rtriGFPl','1GFP','1GFPb','2GFPl1','5GFPa','5GFPa1','5GFPb','5GFPb1','5GFP2la'};
pathn = 'D:/Doc/Sykora/'
for j=1:length(names) 
    %     anti = pt3pro([names{j} '_01.pt3'],[names{j} '_02.pt3'],'antibunch',[],3e9);
    %     [res head] = pt3pro([names{j} '_01.pt3'], [names{j} '_02.pt3'], 'xfcs');
    %     tst = strfind(names{j},'\');
    %     eval(['save ' names{j}(tst(end)+1:end) ' anti res head']);

    load([pathn names{j}])
    
    subplot(2,2,2);
    ind = res.autotime>0; 
    ind(1) = false;
    p = simplex('rigler',[1e-3 1e-4 5e5 1e6],[],[],[],[],res.autotime(ind),mean(res.auto(ind,:),2));
    for k=1:2
        p = simplex('rigler',p,[],[],[],[],res.autotime(ind),mean(res.auto(ind,:),2));
    end
    [err,c] = rigler(p,res.autotime(ind),mean(res.auto(ind,:),2));
    [err,c,z] = rigler(p,c,res.autotime(2:end));
    semilogx(res.autotime(2:end),c(1)+c(2)*z(:,2),'--c',res.autotime(2:end),mean(res.auto(2:end,:),2),'or',res.autotime(2:end),z*c,'b');

    % load([pathn names{j} '_ab'])
    
    subplot(2,2,3);
    pp = simplex('AntiBunch',[0 3 2 8 1],[],[],[],[],anti.autotime,anti.auto,25,1);
    for k=1:2
        pp = simplex('AntiBunch',pp,[],[],[],[],anti.autotime,anti.auto,25,1);
    end
    [err,cc,c0,z,z,z] = AntiBunch(pp,anti.autotime,anti.auto,25,1);
    plot([-anti.autotime(end:-1:2); anti.autotime],[anti.auto(end:-1:2,2); anti.auto(:,1)],'o',[-anti.autotime(end:-1:2); anti.autotime],z(:,1));
    subplot(2,2,4);
    plot([-anti.autotime(end:-1:2); anti.autotime],[anti.auto(end:-1:2,4); anti.auto(:,3)],'o',[-anti.autotime(end:-1:2); anti.autotime],z(:,2));

    h = cc(2)-c0(2);
    delta = cc(3);
    lam = c(2)/sum(c(2:end));
    nemit = lam*h/(lam*h-h+delta);
    npart = c0(2)/h/lam;

    a = sum(cc(1:2))*c(2)/sum(c);
    b = sum(cc(1:2))*(1-c(1)/sum(c))-cc(3);
    nemit = a/(a-b);
    
    subplot(2,2,1)
    cla
    axis off
    text(0.1,0.7,{['\tau_{\itdiff\rm} = ' mnum2str(p(1)*1e3,2,2) ' ms'], ...
        ['\langleemitters\rangle = ' mnum2str(nemit,1,2)], ...
        ['\langleparticles\rangle = ' mnum2str(npart,1,2)]});
    
    subplot(2,2,2)
    xlabel('time [s]');
    ylabel('autocorrelation');
    ax = axis;
    axis([50e-9 3 ax(3) ax(4)])

    subplot(2,2,3)
    ax = axis;
    axis([-112.5 112.5 ax(3) ax(4)])
    xlabel('time [ns]');
    ylabel('autocorrelation');

    subplot(2,2,4)
    axis([-112.5 112.5 ax(3) ax(4)])
    xlabel('time [ns]');
    ylabel('autocorrelation');
    
    eval(['print -dpng ' pathn names{j} '_alt'])
end


