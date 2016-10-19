function [cydyn, hcydyn] = Cy5dynamics(fname,tau,tcspc);

photons = 1e6; % number of photons processed one at a time
smooth = 1; % smoothed autocorrelation
cut = 1; % cut autocorrelatin at 1 mus
close all   

[tmp, tmp, tmp, hcydyn] = tttrRead(fname,[1 1]);
timeunit = hcydyn.GlobClock*1e-9;    
resolution = hcydyn.Resolution;
bin = 1:hcydyn.NChannels;

if nargin==1
    [res0, head] = tttrpro(fname, 'xfcs');
    tau = res0.tau;
    tcspc = res0.tcspc;
    auto0 = res0.auto;
    autotime0 = res0.autotime;
end

t1 = 1:length(tau); t2 = t1; 
t1 = t1(tcspc(:,1)==max(tcspc(:,1)))+25:t1(end-10);
t2 = t2(tcspc(:,2)==max(tcspc(:,2)))+25:t2(end);
t2 = t2(1:length(t1));

p = simplex('expfun',[1 2],[0 0],[],[],[],tau(t1),tcspc(t1,1)); 
taucy = 1./p
p = simplex('expfun',[1 2],[0 0],[],[],[],tau(t2),tcspc(t2,2));
taucy(:,2) = 1./p;

[err, ccy1, zzcy1, zcy(:,1)] = expfun(1./taucy(:,1), tau(t1), tcspc(t1,1)-min(tcspc(t1,1)));
[err, ccy2, zzcy2, zcy(:,2)] = expfun(1./taucy(:,2), tau(t2), tcspc(t2,2)-min(tcspc(t2,2)));

tau1 = tau(t1)'/resolution;
tau2 = tau(t2)'/resolution;
para = [zzcy1(:,1) zzcy2(:,1) zzcy1(:,2) zzcy2(:,2) ones(size(zzcy1,1),2)];

Nsub = 10;
Ncasc = 20;
autotime = cumsum(reshape(repmat(2.^[0 0:Ncasc-1],Nsub,1),(Ncasc+1)*Nsub,1));
autotime = autotime*timeunit;
    
% filter generation
% para = [dye1 channel#1; dye1 channel#2; dye2 channel#1; dye2 channel#2];
for j=1:size(para,2)
    para(:,j) = para(:,j)/sum(para(:,j));
end % normalization
weight = [tcspc(t1,1) tcspc(t2,2)];    
master = (para(:,1:2:end)'*diag(1./weight(:,1))*para(:,1:2:end))\(diag(1./weight(:,1))*para(:,1:2:end))';
master = [master; (para(:,2:2:end)'*diag(1./weight(:,2))*para(:,2:2:end))\(diag(1./weight(:,2))*para(:,2:2:end))']';
    
auto = zeros(length(autotime),size(master,2)/2,size(master,2)/2);
time = 0;
cnt = 0;
tst = 1;
while tst>0
    [tmpy, flv, tmpx, head, tst, overcount, marker] = tttrRead(fname,[cnt+1 photons]);        
    tmpy(marker>0) = [];
    flv(marker>0) = [];
    tmpx(marker>0) = [];        
    if tst>0
        cnt = cnt + tst;
        flv = flv==0;
        ind = (flv & (tmpx<tau1(1) | tau1(end)<tmpx)) | (~flv & (tmpx<tau2(1) | tau2(end)<tmpx));
        tmpy(ind) = [];
        flv(ind) = [];
        tmpx(ind) = [];            
        tmpx(flv) = tmpx(flv) - tau1(1) + 1;
        tmpx(~flv) = tmpx(~flv) - tau2(1) + 1;            
        dt = tmpy(end)-tmpy(1);
        time = time + dt*timeunit;
        flv = double(flv);
        for j=1:size(master,2)/2
            for k=1:size(master,2)/2                
                auto(:,j,k) = auto(:,j,k) + dt*tttr2fcs(tmpy, master(tmpx,j).*flv, ...  
                    master(tmpx,end/2+k).*(1-flv), Ncasc, Nsub, smooth);
            end
        end
        
        subplot('position', [0.925 0.2 0.025 0.6]);
        bar(0,cnt/head.NCounts);
        axis([-0.4 0.4 0 1]);
        set(gca,'xtick',[],'ytick',[]);
        subplot('position', [0.1 0.1 0.7 0.8]);
        semilogx(autotime, auto(:,:,1)/time); 
        drawnow;
    end
end
clf 
if cut
    ind = autotime<5e-7;
    auto(ind,:,:) = [];
    autotime(ind) = [];
end
semilogx(autotime, auto(:,:,1)/time); 
cydyn = struct('autotime', autotime);
cydyn = setfield(cydyn, 'auto', auto/time);
cydyn = setfield(cydyn, 'autotime0', autotime0);
cydyn = setfield(cydyn, 'auto0', auto0);
cydyn = setfield(cydyn, 'tau', head.Resolution*[tau1 tau2]);
cydyn = setfield(cydyn, 'tcspc', weight);
cydyn = setfield(cydyn, 'master', master);    
cydyn = setfield(cydyn, 'weight', weight);
cydyn = setfield(cydyn, 'lifetime', taucy);
cydyn = setfield(cydyn, 'amp', [ccy1 ccy2]);
hcydyn.NCounts = cnt;

