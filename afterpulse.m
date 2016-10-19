% afterpulse

photons = 1e6;
smooth = 1;
timeunit = 1e-7;
name = 'D:\Joerg\Doc\Gregor\040518\NA-B-Atto655_01.t3r';

[res head] = tttrpro(name,'fcs');

tt=10:length(res.tau)-10;

Nsub = 10;
Ncasc = 17;
autotime = cumsum(reshape(repmat(2.^[0 0:Ncasc-1],Nsub,1),(Ncasc+1)*Nsub,1));
autotime = autotime*timeunit;

% filter generation
para = [res.tcspc(tt,1)-min(res.tcspc(tt,1)) res.tcspc(tt,2)-min(res.tcspc(tt,2)) ones(length(tt),2)];
%para = [res.tcspc(tt,:)-ones(length(tt),1)*mean(res.tcspc(end-10:end,:)) ones(length(tt),2)];
for k=1:size(para,2)
    para(:,k) = para(:,k)/sum(para(:,k));
end % normalization
weight = res.tcspc(tt,1:2);    
master = (para(:,1:2:end)'*diag(1./weight(:,1))*para(:,1:2:end))\(diag(1./weight(:,1))*para(:,1:2:end))';
master = [master; (para(:,2:2:end)'*diag(1./weight(:,2))*para(:,2:2:end))\(diag(1./weight(:,2))*para(:,2:2:end))']';
    
tau = res.tau(tt)/head.Resolution;
auto = zeros(length(autotime),4);
crcr = auto;
time = 0;
cnt = 0;
tst = 1;
while tst>0
    [tmpy, flv, tmpx, h, tst] = tttrRead(name,[cnt+1 photons]);        
    if tst>0
        cnt = cnt + tst;
        flv = flv==0;
        ind = tmpx<tau(1) | tau(end)<tmpx;
        tmpy(ind) = [];
        flv(ind) = [];
        tmpx(ind) = [];            
        tmpx = tmpx - tau(1) + 1;
        dt = tmpy(end)-tmpy(1);
        time = time + dt*timeunit;
        flv = double(flv);
        auto(:,1) = auto(:,1) + dt*tttr2fcs(tmpy, master(tmpx,1).*flv, master(tmpx,1).*flv, Ncasc, Nsub, smooth);
        auto(:,2) = auto(:,2) + dt*tttr2fcs(tmpy, master(tmpx,3).*(1-flv), master(tmpx,3).*(1-flv), Ncasc, Nsub, smooth);        
        auto(:,3) = auto(:,3) + dt*tttr2fcs(tmpy, master(tmpx,2).*flv, master(tmpx,2).*flv, Ncasc, Nsub, smooth);
        auto(:,4) = auto(:,4) + dt*tttr2fcs(tmpy, master(tmpx,4).*(1-flv), master(tmpx,4).*(1-flv), Ncasc, Nsub, smooth);        
        crcr(:,1) = crcr(:,1) + dt*tttr2fcs(tmpy, master(tmpx,1).*flv, master(tmpx,2).*flv, Ncasc, Nsub, smooth);
        crcr(:,2) = crcr(:,2) + dt*tttr2fcs(tmpy, master(tmpx,3).*(1-flv), master(tmpx,4).*(1-flv), Ncasc, Nsub, smooth);        
        crcr(:,3) = crcr(:,3) + dt*tttr2fcs(tmpy, master(tmpx,2).*flv, master(tmpx,1).*flv, Ncasc, Nsub, smooth);
        crcr(:,4) = crcr(:,4) + dt*tttr2fcs(tmpy, master(tmpx,4).*(1-flv), master(tmpx,3).*(1-flv), Ncasc, Nsub, smooth);        
        
        subplot('position', [0.925 0.2 0.025 0.6]);
        bar(0,cnt/h.NCounts);
        axis([-0.4 0.4 0 1]);
        set(gca,'xtick',[],'ytick',[]);
        subplot('position', [0.1 0.1 0.7 0.8]);
        semilogx(autotime, auto/time); 
        drawnow;
    end
end
auto = auto/time;
crcr = crcr/time;

break

load afterpulse

tmp=auto(10:end,2)/auto(end,2)*res.auto(end,1,1);
semilogx(res.autotime,res.auto(:,1,1)/res.auto(end,1,1),'ob',res.autotime,tmp/tmp(end),'r')
axis tight
xlabel('time [s]');
ylabel('autocorrelation')

semilogy(tau,weight(:,2),'r');
xlabel('time [ns]');
ylabel('photon counts')

semilogy(tau,para(:,2),'r',tau,para(:,4),'b');
xlabel('time [ns]');
ylabel('master curves [a.u.]')

