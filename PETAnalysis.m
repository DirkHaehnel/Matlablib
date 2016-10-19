function PETAnalysis(dname, attoname, petname)

close all

load([dname attoname])
[y,t] = FCSCrossReadCW(res);
ind = t>1e-8 & t<1;

pdye = Simplex('Rigler',[0.004 5 2],[],[],[],[],t(ind),abs(mean(sum(y(ind,:,:),3),2)),2,1);
for j=1:10
    pdye = Simplex('Rigler',pdye,[],[],[],[],t(ind),abs(mean(sum(y(ind,:,:),3),2)),2,1);
end
[~, cdye, zdye, zzdye] = Rigler(pdye,t(ind),abs(mean(sum(y(ind,:,:),3),2)),2,1);

fnames = dir([dname petname]);
ind1 = t>1e-8 & t<1e0;
for j=1:length(fnames)
    close
    
    load([dname fnames(j).name])
    
    savename = strrep(fnames(j).name(1:end-4),'.','-');
    
    [y,t] = FCSCrossReadCW(res);
    
    p1 = Simplex('AffineFit',[1 1e-8], [], [], [], [], t(ind1), abs(mean(mean(y(ind1,:,:),3),2)), t(ind), zdye/max(zdye), [], 1, 1);
    for k=1:10
       p1 = Simplex('AffineFit',p1, [], [], [], [], t(ind1), abs(mean(mean(y(ind1,:,:),3),2)), t(ind), zdye/max(zdye), [], 1, 1); 
    end
    [~, c1, z1, zz1, back1] = AffineFit(p1, t(ind1), abs(mean(mean(y(ind1,:,:),3),2)), t(ind), zdye/max(zdye), [], 1, 1);
    tau_on = 1/(sum(c1(2))/sum(c1(2:end))/p1(end))
    tau_off = 1/(sum(c1(3:end))/sum(c1(2:end))/p1(end))
    
    subplot(4,1,1:3)
    semilogx(t(ind1),abs(mean(mean(y(ind1,:,:),3),2)),'o',t(ind1),z1)
    ylabel('autocorrelaton')
    title(fnames(j).name(1:end-4),'Interpreter','none')
    ax = axis;
    text(ax(1)+1e-3*diff(ax(1:2)),ax(3)+0.7*diff(ax(3:4)),{['\tau_{on} = ' int2str(round(tau_on*1e9)) ' ns'],...
        ['\tau_{off} = ' int2str(round(tau_off*1e9)) ' ns'],[],...
        ['\tau_{total} = ' int2str(round(p1(end)*1e9)) ' ns'],[],...
        ['\tau_{\it D} = ' mnum2str(p1(1),1,1) ' \tau_{Atto}']})
    subplot(4,1,4)
    semilogx(t(ind1),abs(mean(mean(y(ind1,:,:),3),2))-z1)
    ylabel('residuals')
    xlabel('lag time (s)')
    
    %eval(['print -r300 -dpng ' '''' dname fnames(j).name(1:end-4) '''' 'A'])
    eval(['print -r300 -dpng ' '''' dname savename ''''])
    
    close
    
    plot(cumsum(res.time(sum(res.rate,2)>0)),res.rate(sum(res.rate,2)>0,:))
    xlabel('time (s)');
    ylabel('count rate');
    title(fnames(j).name(1:end-4),'Interpreter','none')

    eval(['print -r300 -dpng ' '''' dname savename '''' '_Rate'])
    
end

