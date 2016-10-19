function res = BurstAnalysis(fname,a)

close all
photons = 1e7;
wnd = 15; % width of time window üper pulse in ns

bld = true;
flag = true;
cnt = 0;
lastind = 0;
res.bs = [];
res.sync = [];
res.tcspc = [];
res.chan = [];
res.ind = [];

while flag
    [sync, tcspc, chan, ~, num, ~, head] = HT3_Read(fname,[cnt photons]); % Im Argument von HT3_Read cnt+photons zu photons geändert
    if num==0
        flag = false;
    else
        cnt = cnt + num;
        if bld
            % extract and align lifetime curves; excitation blue/blue/red/red;
            % detection red/red/blue/blue
            tautime = min(tcspc):max(tcspc);
            lt = zeros(numel(tautime),4);
            for j=1:4
                lt(:,j)=mHist(tcspc(chan==j-1),tautime);
            end
            semilogy(head.Resolution*tautime(1:end-10),lt(1:end-10,:))
            xlabel('time (ns)');
            ylabel('fluorescence (cnts)');
            legend({'red 1', 'red 2', 'blue 1', 'blue 2'})
            eval(['print -dpng -r300 ''' fname(1:end-4) '_LifetimeRaw'''])
            res.head = head;
            
            if nargin==1 || isempty(a)
                a = zeros(5,2);
                pm = zeros(2,4);
                semilogy(tautime(1:end-10),lt(1:end-10,1))
                a(1:2,:) = ginput(2);
                for j=2:4
                    semilogy(tautime(1:end-10),lt(1:end-10,j))
                    a(1+j,:) = ginput(1);
                end
                a = round(a(:,1));
            end
            res.a = a';
            tauwin = 1:round(wnd/head.Resolution);
            lw = zeros(numel(tauwin),2,4);
            for k=1:2
                pm(k,1) = a(k);
            end
            for j=2:4
                pm(:,j) = a(1+j) + pm(:,1) - pm(1,1);
            end
            for j=1:4
                for k=1:2
                    lw(:,k,j) = lt(pm(k,j)+tauwin,j);
                end
            end
            tw = tauwin;
            res.tau = tw;
            res.pm = pm;
            
            close
            subplot(211)
            semilogy(tw,lw(:,1,1),'b',tw,lw(:,1,2),'c:',tw,lw(:,2,1),'r',tw,lw(:,2,2),'m:');
            xlabel('time (ns)');
            ylabel('red fluorescence');
            axis tight
            legend({'BxRd1', 'BxRd2', 'RxRd1', 'RxRd2'},-1)
            subplot(212)
            semilogy(tw,lw(:,1,3),'b',tw,lw(:,1,4),'c:',tw,lw(:,2,3),'r',tw,lw(:,2,4),'m:');
            xlabel('time (ns)');
            ylabel('blue fluorescence');
            axis tight
            legend({'BxBd1', 'BxBd2', 'RxBd1', 'RxBd2'},-1)
            eval(['print -dpng -r300 ''' fname(1:end-4) '_Lifetime'''])
            bld = false;
        end
        
        % determine excitation and detection channels
        RxRd1 = pm(2,1)<tcspc & tcspc<=(pm(2,1)+tauwin(end)) & chan==0;
        RxRd2 = pm(2,2)<tcspc & tcspc<=(pm(2,2)+tauwin(end)) & chan==1; 
        RxBd1 = pm(2,3)<tcspc & tcspc<=(pm(2,3)+tauwin(end)) & chan==2; 
        RxBd2 = pm(2,4)<tcspc & tcspc<=(pm(2,4)+tauwin(end)) & chan==3; 
        BxRd1 = pm(1,1)<tcspc & tcspc<=(pm(1,1)+tauwin(end)) & chan==0; 
        BxRd2 = pm(1,2)<tcspc & tcspc<=(pm(1,2)+tauwin(end)) & chan==1; 
        BxBd1 = pm(1,3)<tcspc & tcspc<=(pm(1,3)+tauwin(end)) & chan==2; 
        BxBd2 = pm(1,4)<tcspc & tcspc<=(pm(1,4)+tauwin(end)) & chan==3;
        tcspc(RxRd1) = tcspc(RxRd1) - pm(2,1);
        tcspc(RxRd2) = tcspc(RxRd2) - pm(2,2);
        tcspc(RxBd1) = tcspc(RxBd1) - pm(2,3);
        tcspc(RxBd2) = tcspc(RxBd2) - pm(2,4);
        tcspc(BxRd1) = tcspc(BxRd1) - pm(1,1);
        tcspc(BxRd2) = tcspc(BxRd2) - pm(1,2);
        tcspc(BxBd1) = tcspc(BxBd1) - pm(1,3);
        tcspc(BxBd2) = tcspc(BxBd2) - pm(1,4);
        
        % filter all photons which lie within chosen time windows
        chan(BxRd1 | BxRd2 | BxBd1 | BxBd2) = chan(BxRd1 | BxRd2 | BxBd1 | BxBd2) + 4;
        ind = RxRd1 | RxRd2 | RxBd1 | RxBd2 | BxRd1 | BxRd2 | BxBd1 | BxBd2;
        sync = sync(ind);
        tcspc = tcspc(ind);
        chan = chan(ind);
        
        % binning
        % order: [RxRd1 RxRd2 RxBd1 RxBd2 BxRd1 BxRd2 BxBd1 BxBd2]
        bb = tttr2bin(sync, round(head.SyncRate*1e-4), chan);
        
        % FRET burst search
        minbs = 10; % minimum burst size
        m = 10; % half-width smoothing window
        cap = 0;
        mm = 100; % background smoothing window
        fac = 1.5;
        [bs, t1, t2] = bin2burst(bb, minbs, m, cap, mm, fac);
        if ~isempty(bs)
            res.bs = [res.bs; bs];
            
            fretind = sum(bs(:,5:6),2)>10; % bursts with more than 10 photons upon BxRd
            
            % show accumulated FRET histogram
            clf
            mHist(sum(res.bs(fretind,5:6),2)./(sum(res.bs(fretind,5:6),2)+sum(res.bs(fretind,7:8),2)),0:0.05:1)
            xlabel('FRET efficiency'), ylabel('frequency');
            drawnow
            
            tt = cumsum(sum(bb,2));
            ind1 = zeros(numel(sync),1);
            ind2 = ind1;
            ind1(tt(t1-1)+1) = 1;
            ind2(tt(t2)+1) = 1;
            ind1 = (cumsum(ind1)-cumsum(ind2)>0).*cumsum(ind1);
            res.sync = [res.sync; sync(ind1>0)];
            res.tcspc = [res.tcspc; tcspc(ind1>0)];
            res.chan = [res.chan; chan(ind1>0)];
            res.ind = [res.ind; lastind + ind1(ind1>0)];
            lastind = res.ind(end);
        end
    end
    
    tt = 1:numel(res.ind);
    res.t1 = tt(diff([0; res.ind(1:end)])==1);
    res.t2 = tt(diff([res.ind; max(res.ind)+1])==1);
    
    eval(['save ''' fname(1:end-4) '_results''' ' res'])
end

return


mHist2(sum(bs(:,5:6),2),sum(bs(:,7:8),2),1:1:100,1:1:150); xlabel('FRET'); ylabel('Donor'); shading flat

mHist(sum(bs(fretind,5:6),2)./(sum(bs(fretind,5:6),2)+sum(bs(fretind,7:8),2)),0:0.05:1)
xlabel('FRET efficiency'), ylabel('frequency');

close all; 
ind = res.tw<10;
close all; tau_acceptor_fret=Simplex('ExpFun',[0.5 1],[],[],[],[],res.tw(ind),res.acceptor_fret(ind,1),2)
xlabel('time (ns)'), ylabel('acceptor');
axis tight
ax = axis;
text(ax(1)+0.8*diff(ax(1:2)),ax(3)+0.5*diff(ax(3:4)),{['\tau_1 = ' mnum2str(tau_acceptor_fret(1),1,2)],['\tau_2 = ' mnum2str(tau_acceptor_fret(2),1,2)]})
eval(['print -dpng -r300 ''' fname(1:end-4) '_AcceptorFRETLifetime'''])

close all; 
ind = 0.1<res.tw & res.tw<10;
tau_donor_fret=Simplex('ExpFun',[0.1 0.5 3],[],[],[],[],res.tw(ind),res.donor_fret(ind,1),2)
xlabel('time (ns)'), ylabel('donor');
axis tight
ax = axis;
text(ax(1)+0.8*diff(ax(1:2)),ax(3)+0.5*diff(ax(3:4)),{['\tau_1 = ' mnum2str(tau_donor_fret(1),1,2)],['\tau_2 = ' mnum2str(tau_donor_fret(2),1,2)],['\tau_3 = ' mnum2str(tau_donor_fret(3),1,2)]})
eval(['print -dpng -r300 ''' fname(1:end-4) '_DonorFRETLifetime'''])

close all; 
ind = 0.1<res.tw & res.tw<10;
tau_donor_only=Simplex('ExpFun',[0.1 3],[],[],[],[],res.tw(ind),res.donor_only(ind,1),2)
xlabel('time (ns)'), ylabel('donor only');
axis tight
ax = axis;
text(ax(1)+0.8*diff(ax(1:2)),ax(3)+0.5*diff(ax(3:4)),{['\tau_1 = ' mnum2str(tau_donor_only(1),1,2)],['\tau_2 = ' mnum2str(tau_donor_only(2),1,2)]})
eval(['print -dpng -r300 ''' fname(1:end-4) '_DonorOnlyLifetime'''])

close all; 
ind = res.tw<10;
close all; tau_acceptor=Simplex('ExpFun',[0.1 1],[],[],[],[],res.tw(ind),res.acceptor(ind,1),2)
xlabel('time (ns)'), ylabel('acceptor');
axis tight
ax = axis;
text(ax(1)+0.8*diff(ax(1:2)),ax(3)+0.5*diff(ax(3:4)),{['\tau_1 = ' mnum2str(tau_acceptor(1),1,2)],['\tau_2 = ' mnum2str(tau_acceptor(2),1,2)]})
eval(['print -dpng -r300 ''' fname(1:end-4) '_AcceptorDirectLifetime'''])

return

fname = 'm:\MTusers\Qui\120717_HTH-EMBL\Point_005.ht3';
res=BurstAnalysis(fname,[1354 13950 1296 989 1046]);