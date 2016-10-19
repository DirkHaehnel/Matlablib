function [res, head] = PicosecondCorrelation(name)

automaxtime = 3e-6;
photons = 1e6; 
close all

if strcmp(name(end-2:end),'ht3')
	head = ht3v2read(name);
    if isempty(head)
        if nargin<3 || isempty(para)
            Timeunit = 50e-9;
            Resolution = 2e-12;
        else
            Timeunit = para(1)*1e-9;
            Resolution = para(2)*1e-9;
        end
        NChannels  = ceil(Timeunit/Resolution);
        tmp = dir(name);
        NCounts = tmp.bytes/4;
        dind = [0 1];
    else
        Timeunit = 1/(head.SyncRate);
        Resolution = head.Resolution*1e-12;
        NChannels  = ceil(Timeunit/Resolution);
        tmp = dir(name);
        NCounts = tmp.bytes/4;
        dind = [0 1];
    end
else
    disp('Not a valid file!');
    return
end

autotime = (0:automaxtime/Resolution)*Resolution;
auto = zeros(length(autotime),2);

clf
tcspc = zeros(NChannels,2);
rate = [];
time = 0;
jj = 1;
cnt = 0;
tst = 1;
while tst>0
    [sync, micro, chan, bla, tst] = ht3read(name, [cnt+1, photons]);
    if tst>0
        cnt = cnt + tst;

        tcspc(:,1) = tcspc(:,1) + mHist(micro(chan==dind(1)),1:NChannels);
        tcspc(:,2) = tcspc(:,2) + mHist(micro(chan==dind(2)),1:NChannels);

        tim = sync*Timeunit + micro*Resolution;
        time(jj) = tim(end) - tim(1);
        
        for k=1:ceil(automaxtime/150e-9)
            tmp = tim(k+1:end)-tim(1:end-k);
            ind = chan(1:end-k)==dind(1) & chan(k+1:end)==dind(2);
            if sum(ind)>1
                auto(:,1) = auto(:,1) + mHist(tmp(ind),autotime);
            end
            ind = chan(1:end-k)==dind(2) & chan(k+1:end)==dind(1);
            if sum(ind)>1
                auto(:,2) = auto(:,2) + mHist(tmp(ind),autotime);
            end
        end

        rate(jj,1) = sum(chan==dind(1))/time(jj);
        rate(jj,2) = sum(chan==dind(2))/time(jj);

        subplot('position', [0.925 0.2 0.025 0.6]);
        bar(0,cnt/NCounts);
        axis([-0.4 0.4 0 1]);
        set(gca,'xtick',[],'ytick',[]);
        subplot('position', [0.1 0.1 0.7 0.8]);
        plot(-autotime, auto(:,2), autotime, auto(:,1));
        drawnow;
    end
    jj = jj+1;
end

clf
plot(-autotime, auto(:,2), autotime, auto(:,1));

res.tau = (1:NChannels)*Resolution;
res.tcspc = tcspc;
res.autotime = autotime;
res.auto = auto;
res.rate = rate;
res.time = time;


return

width = 25e3;
t = mean(reshape(res.autotime(1:floor(length(res.autotime)/width)*width),width,floor(length(res.autotime)/width)));
clear a
a(:,1) = sum(reshape(res.auto(1:floor(length(res.autotime)/width)*width,1),width,floor(length(res.autotime)/width)))';
a(:,2) = sum(reshape(res.auto(1:floor(length(res.autotime)/width)*width,2),width,floor(length(res.autotime)/width)))';
plot(-1e9*t,a(:,2),1e9*t,a(:,1))
xlabel('time (ns)'); ylabel('correlation')