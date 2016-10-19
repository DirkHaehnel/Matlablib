% AntibunchAnalysis

pathn = 'D:\Joerg\Doc\Fcs\Antibunching\';
names = {'Atto655'}; 
bild = false;
for j=1:length(names) 
    anti = pt3pro([pathn names{j} '_01.pt3'],[pathn names{j} '_02.pt3'],'antibunchm',[],3e9);
    [res head] = pt3pro([pathn names{j} '_01.pt3'], [pathn names{j} '_02.pt3'], 'xfcsm');
    eval(['save ' pathn names{j} ' anti res head']);

    return 
    
    load([pathn names{j}])
    
    if bild
        subplot(2,2,2);
        %ind = res.autotime<5e-6 | res.autotime>1;
        ind = res.autotime>0;
        ind(1) = false;
        p = simplex('rigler',[1e-3 1e-4 1e5],[],[],[],[],res.autotime(ind),mean(res.auto(ind,:),2));
        for k=1:2
            p = simplex('rigler',p,[],[],[],[],res.autotime(ind),mean(res.auto(ind,:),2));
        end
        [err,c] = rigler(p,res.autotime(ind),mean(res.auto(ind,:),2));
        [err,c,z] = rigler(p,c,res.autotime(2:end));
        semilogx(res.autotime(2:end),c(1)+c(2)*z(:,2),'--c',res.autotime(2:end),mean(res.auto(2:end,:),2),'or',res.autotime(2:end),z*c,'b');

        subplot(2,2,3);
        pp = simplex('AntiBunch',[0 3 2 8 1],[],[],[],[],anti.autotime,anti.auto,25,1);
        for k=1:2
            pp = simplex('AntiBunch',pp,[],[],[],[],anti.autotime,anti.auto,25,1);
        end
        [err,cc,c0,z,z,z] = AntiBunch(pp,anti.autotime,anti.auto,25,1);
        %     pp = simplex('AntiBunchFun',[0 2 8 1],[],[],[],[],anti.autotime,anti.auto,ctau,tau,25,1);
        %     for k=1:2
        %         pp = simplex('AntiBunchFun',pp,[],[],[],[],anti.autotime,anti.auto,ctau,tau,25,1);
        %     end
        % [err,cc,c0,z,z,z] = AntiBunchFun(pp,anti.autotime,anti.auto,ctau,tau,25,1);
        plot([-anti.autotime(end:-1:2); anti.autotime],[anti.auto(end:-1:2,2); anti.auto(:,1)],'o',[-anti.autotime(end:-1:2); anti.autotime],z(:,1));
        subplot(2,2,4);
        plot([-anti.autotime(end:-1:2); anti.autotime],[anti.auto(end:-1:2,4); anti.auto(:,3)],'o',[-anti.autotime(end:-1:2); anti.autotime],z(:,2));

        h = cc(2)-c0(2);
        delta = cc(3);
        lam = c(2)/sum(c(2:end));
        nemit = lam*h/(lam*h-h+delta);
        npart = c0(2)/h/lam;

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
        axis([50e-9 1 ax(3) ax(4)])

        subplot(2,2,3)
        ax = axis;
        axis([-112.5 112.5 ax(3) ax(4)])
        xlabel('time [ns]');
        ylabel('autocorrelation');

        subplot(2,2,4)
        axis([-112.5 112.5 ax(3) ax(4)])
        xlabel('time [ns]');
        ylabel('autocorrelation');
    else
        close all
        %ind = res.autotime<5e-6 | res.autotime>1;
        ind = res.autotime>0;
        ind(1) = false;
        p = simplex('rigler',[1e-2 1e-3 1e7 1e6],[],[],[],[],res.autotime(ind),mean(res.auto(ind,:),2));
        for k=1:5
            p = simplex('rigler',p,[],[],[],[],res.autotime(ind),mean(res.auto(ind,:),2));
        end
        [err,c] = rigler(p,res.autotime(ind),mean(res.auto(ind,:),2));
        [err,c,z] = rigler(p,c,res.autotime(2:end));
        semilogx(res.autotime(2:end),c(1)+c(2)*z(:,2),'--c',res.autotime(2:end),mean(res.auto(2:end,:),2),'or',res.autotime(2:end),z*c,'b');
        xlabel('time [s]');
        ylabel('autocorrelation');
        ax = axis;
        axis([50e-9 1 ax(3) ax(4)])

        figure
        pp = simplex('AntiBunch',[0 3 2 8 1],[],[],[],[],anti.autotime,anti.auto,25,1);
        for k=1:2
            pp = simplex('AntiBunch',pp,[],[],[],[],anti.autotime,anti.auto,25,1);
        end
        [err,cc,c0,z,z,z] = AntiBunch(pp,anti.autotime,anti.auto,25,1);
        %     pp = simplex('AntiBunchFun',[0 2 8 1],[],[],[],[],anti.autotime,anti.auto,ctau,tau,25,1);
        %     for k=1:2
        %         pp = simplex('AntiBunchFun',pp,[],[],[],[],anti.autotime,anti.auto,ctau,tau,25,1);
        %     end
        % [err,cc,c0,z,z,z] = AntiBunchFun(pp,anti.autotime,anti.auto,ctau,tau,25,1);
        plot([-anti.autotime(end:-1:2); anti.autotime],[anti.auto(end:-1:2,2); anti.auto(:,1)],'o',[-anti.autotime(end:-1:2); anti.autotime],z(:,1));
        ax = axis;
        axis([-112.5 112.5 ax(3) ax(4)])
        xlabel('time [ns]');
        ylabel('autocorrelation');
        figure
        plot([-anti.autotime(end:-1:2); anti.autotime],[anti.auto(end:-1:2,4); anti.auto(:,3)],'o',[-anti.autotime(end:-1:2); anti.autotime],z(:,2));
        axis([-112.5 112.5 ax(3) ax(4)])
        xlabel('time [ns]');
        ylabel('autocorrelation');

        h = cc(2)-c0(2);
        delta = cc(3);
        lam = c(2)/sum(c(2:end));
        nemit = lam*h/(lam*h-h+delta);
        npart = c0(2)/h/lam;
    end
    % eval(['print -dpng ' names{j}(tst(end)+1:end)])
    % eval(['print -r300 -dpng ' pathn names{j}])
end

return

for jf=1:length(files)

    names = files{jf};

    subplot(211)

    clear res flcs
    load(names{1})
    if ~exist('flcs')==1
        p = simplex('rigler',[1e-2 1e-3 1e6],[],[],[],[],res.autotime(2:end),mean(res.auto(2:end,:),2));
        for k=1:2
            p = simplex('rigler',p,[],[],[],[],res.autotime(2:end),mean(res.auto(2:end,:),2));
        end
        [err,c] = rigler(p,res.autotime(2:end),mean(res.auto(2:end,:),2),[],1);
    else
        p = simplex('rigler',[1e-2 1e-3 1e6],[],[],[],[],flcs.autotime(2:end),flcs.auto(2:end,1,1));
        for k=1:2
            p = simplex('rigler',p,[],[],[],[],res.autotime(2:end),flcs.auto(2:end,1,1));
        end
        [err,c] = rigler(p,res.autotime(2:end),flcs.auto(2:end,1,1),[],1);
    end
    xlabel('time [s]');
    ylabel('autocorrelation');
    ax = axis;

    subplot(212)

    load(names{2})
    pp = simplex('AntiBunch',[0 3 2 8 1],[],[],[],[],res.autotime,res.auto(:,[1 2]),25,1);
    for k=1:2
        pp = simplex('AntiBunch',pp,[],[],[],[],res.autotime,res.auto(:,[1 2]),25,1);
    end
    [err,nemit,sig,cc] = AntiBunch(pp,res.autotime,res.auto(:,[1 2]),25,1,1);

    ind = res.tau>3.5 & res.tau<24;
    tau1 = simplex('expfun',3,[],[],[],[],res.tau(ind),res.tcspc(ind,1));
    tau2 = simplex('expfun',3,[],[],[],[],res.tau(ind),res.tcspc(ind,2));
    [err, ctau1] = expfun(tau1,res.tau(ind),res.tcspc(ind,1));
    [err, ctau2] = expfun(tau2,res.tau(ind),res.tcspc(ind,2));
    sig = sum(ctau1(2:end))/sum(ctau1)*sum(ctau2(2:end))/sum(ctau2);

    npart = c(1)/c(2)

    nemit = cc(2)/cc(3)*c(2)/sum(c)
    
    ax = axis;
    axis([-112.5 112.5 ax(3) ax(4)])
    text(0,ax(3)+0.8*(ax(4)-ax(3)),{['\langleemitters\rangle = ' mnum2str(nemit(1),1,2)]})
    xlabel('time [ns]');
    ylabel('autocorrelation');

    subplot(211)
    ax = axis;
    axis([50e-9 3 ax(3) ax(4)])
    text(0.01,ax(3)+0.8*(ax(4)-ax(3)),{['\langleparticles\rangle = ' mnum2str(npart(1),1,2)],...
        ['\tau_{\itdiff\rm} = ' mnum2str(p(1)*1e3,2,2) ' ms']})

    title(names{1},'interpreter','none')

    eval(['print -dpng ' names{1} '_a'])

    subplot(212)
    [err,nemit,sig] = AntiBunch(pp,res.autotime,res.auto(:,[1 2]),25,1);
    npart = (c(1)-sig*sum(c))/sum(c(2:end))

    ax = axis;
    axis([-112.5 112.5 ax(3) ax(4)])
    text(0,ax(3)+0.8*(ax(4)-ax(3)),{['\langleemitters\rangle = ' mnum2str(nemit(1),1,2)],...
        ['\langleemitters per particle\rangle = ' mnum2str(nemit(1)/npart(1),1,2)]})
    xlabel('time [ns]');
    ylabel('autocorrelation');

%     clf 
%     semilogy(res.tau,res.tcspc);
%     ax = axis;
%     text(15,ax(3)+0.1*(ax(4)-ax(3)),{['\tau_1 = ' mnum2str(tau1,1,2) ' ns'],...
%         ['\tau_2 = ' mnum2str(tau2,1,2) ' ns']})
%     xlabel('time [ns]');
%     ylabel('photon counts');
%     title([names{1} ' lifetime'],'interpreter','none')
%     eval(['print -dpng ' names{1} '_lifetime'])

end

return

load(names{2})
pp = simplex('AntiBunch',[0 1 2 8 4],[],[],[],[],res.autotime,res.auto,25,1); 
for k=1:2
    pp = simplex('AntiBunch',pp,[],[],[],[],res.autotime,res.auto,25); 
end
[err,nemit,sig] = AntiBunch(pp,res.autotime,res.auto,25); 

