function [res, head] = SingleFocusCW2FCS(name, mintime, maxtime, Nsub, photons)

if nargin<3 || isempty(maxtime)
    maxtime = 3; % maxtime in s
end

if nargin<4 || isempty(Nsub)
    Nsub = 10; % subdivision of lag times
end

if nargin<5 || isempty(photons)
    photons = 1e6; % number of photons processed one at a time
end

head = ht3v2read(name);

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
    
microresolution = head.Resolution;          % resolution in ns
syncperiod = (1 / head.SyncRate) * 1e9;     % syncperiod in ns

if nargin<2 || isempty(mintime)
    mintime = head.Resolution*1e-9; % in s
end

Ncasc = ceil(log2(maxtime/mintime/Nsub));    
autotime = cumsum(reshape(repmat(2.^(0:Ncasc-1),Nsub,1),Ncasc*Nsub,1));
autotime = autotime*mintime;      % in s

NCounts = head.Records;

[~, ~, chan, special] = ht3v2read(name, [1, photons]);
tchan = chan(special == 0);
dind = unique(tchan(1:1e4));
ndet = length(dind);
rev_dind = zeros(8,1);
for j=1:length(dind)
    rev_dind(dind(j)+1) = j;
end

auto = 0*autotime;
auto = repmat(auto,[1 ndet ndet ceil(NCounts/photons)]);

clf;
tcspc = zeros(2^15, ndet);
rate = zeros(ceil(NCounts/photons), ndet);
time = zeros(ceil(NCounts/photons), 1);
jj = 1;
cnt = 0;
tst = 1;

while tst>0
    [sync, micro, chan, special, tst] = ht3v2read(name, [cnt+1, photons]);
    cnt = cnt + tst;
    if (tst>0.9*photons) && (~isempty(sync))
        ind = special==1 | ~ismember(chan,dind);
        sync(ind) = [];
        micro(ind) = [];
        chan(ind) = [];
        
        for j=1:length(chan)
            dd = rev_dind(chan(j)+1);
            tcspc(micro(j)+1, dd) = tcspc(micro(j)+1, dd) + 1;
        end

        tim = 1e-9*(sync*syncperiod + micro*microresolution);  % photon arrival time in ns
        time(jj) = (tim(end) - tim(1));            % length in s
        
%         for det1 = 1:ndet
%             rate(jj,det1) = sum(chan==dind(det1))/time(jj); % rate per second
%         end
%         rate(jj,:) = mean(tttr2bin(tim,autotime(end),chan))/autotime(end);
%         [sum(chan==0) sum(chan==1)]/time(jj)./(mean(tttr2bin(tim,autotime(end),chan))/autotime(end))
        
        flvp = zeros(length(chan),ndet);
        for p=1:ndet
            flvp(:,p) = chan==dind(p);
        end
        rate(jj,:) = sum(flvp)/time(jj); % rate per second
        auto(:,:,:,jj) = tttr2xfcs(round(tim/mintime), flvp, Ncasc, Nsub) / mintime;

        for det1 = 1:ndet
            for det2 = 1:ndet
                auto(:,det1,det2,jj) = ((auto(:,det1,det2,jj) / time(jj) / (rate(jj,det1) * rate(jj,det2)))) - 1;
            end
        end
        
        subplot('position', [0.925 0.2 0.025 0.6]);
        bar(0,cnt/NCounts);
        axis([-0.4 0.4 0 1]);
        set(gca,'xtick',[],'ytick',[]);
        subplot('position', [0.1 0.1 0.7 0.8]);
        
        [y, t] = FCSCrossReadCW(struct('auto',auto,'autotime',autotime));
        semilogx(t, mean(y(:,:,1:jj),3));
        %ylim([0 max(max(mean(y(:,:,1:jj), 3)))*1.05]);
        drawnow;
    end
    jj = jj+1;
    save test autotime auto
end

clf

[y, t] = FCSCrossReadCW(struct('auto',auto,'autotime',autotime));
semilogx(t, mean(y,3));

tau = (0:2^15-1)*microresolution;
tau = tau';
tau(any(tcspc==0,2),:)=[];
tcspc(any(tcspc==0,2),:)=[];

res.tau = tau;
res.tcspc = tcspc;

res.autotime = autotime;
res.auto = auto;

res.rate = rate;
res.time = time;



