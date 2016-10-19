function [res, head] = SingleFocus2Triplet(name, timebin);

photons = 1e6; % number of photons processed one at a time

close all   

if nargin<2 | isempty(timebin)
    timebin=100; 
end

[tmp, tmp, tmp, head] = tttrRead(name,[1 1]);
timeunit = head.GlobClock*1e-9;    

autotime = (1:timebin);
autotime = autotime*timeunit;
order = floor(timebin/10);
auto = zeros(length(autotime),order,2);

% xtmp = zeros(head.NChannels,2);
% tst = 1;
% cnt = 0;
% bin = 1:head.NChannels;
% while tst>0
%     [tmpy, flv, tmpx, head, tst] = tttrRead(name,[cnt+1 photons]);  
%     cnt = cnt+tst;
%     if tst>0
%         xtmp(:,1) = xtmp(:,1) + mhist(tmpx(flv==0),bin);
%         xtmp(:,2) = xtmp(:,2) + mhist(tmpx(flv==1),bin);
%     end
% end
% ind = sum(xtmp,2)==0;
% bin(ind) = [];
% xtmp(ind,:) = [];
% bin([1:30 end-30:end])=[];
% xtmp([1:30 end-30:end],:) = [];
% semilogy(bin,xtmp); drawnow;
% 
% bin = bin';
% win = (1:round(length(bin)))-30;
% 
% tcspc = zeros(length(win),2); tau = tcspc;
% [pos pos] = max(xtmp(:,1));
% ind = pos + win;
% if ind(1)<=0 
%     win(1:sum(ind<=0)) = [];
%     tcspc(1:sum(ind<=0),:) = [];
%     tau(1:sum(ind<=0),:) = [];        
%     ind = pos + win;
% end
% if ind(end)>length(xtmp(:,1)) 
%     win(end-sum(ind>length(xtmp(:,1))):end) = [];
%     tcspc(end-sum(ind>length(xtmp(:,1))):end,:) = [];
%     tau(end-sum(ind>length(xtmp(:,1))):end,:) = [];
%     ind = pos + win;
% end
% tcspc(:,1) = xtmp(ind,1);
% tau(:,1) = bin(ind);
% [pos pos] = max(xtmp(:,2));
% ind = pos + win;
% if ind(1)<=0 
%     win(1:sum(ind<=0)) = [];
%     tcspc(1:sum(ind<=0),:) = [];
%     tau(1:sum(ind<=0),:) = [];
%     ind = pos + win;
% end
% if ind(end)>length(xtmp(:,2)) 
%     win(end-sum(ind>length(xtmp(:,2))):end) = [];
%     tcspc(end-sum(ind>length(xtmp(:,2))):end,:) = [];
%     tau(end-sum(ind>length(xtmp(:,2))):end,:) = [];
%     ind = pos + win;
% end
% tcspc(:,2) = xtmp(ind,2);
% tau(:,2) = bin(ind);
% 
% tcspc = tcspc-min(min(tcspc));
% semilogy(tau,tcspc); drawnow
% res = struct('tau', head.Resolution*tau);
% res = setfield(res, 'tcspc', tcspc);

time = 0;
cnt = 0;
tst = 1;
t = 1:size(auto,1);
while tst>0
    [tmpy, flv, tmpx, head, tst] = tttrRead(name,[cnt+1 photons]);        
    if tst>0
        cnt = cnt + tst
        dt = tmpy(end)-tmpy(1);
        time = time + dt;
        for j=1:length(tmpy)-order
            tmp = tmpy(j+1:j+order)-tmpy(j);
            flag = flv(j+1:j+order);
            s = tmp<timebin;
            if flv(j)==0
                ind = s & flag==1;
                if sum(ind)>0
                    auto(tmp(ind)+1,sum(s),1) = auto(tmp(ind)+1,sum(s),1) + 1;
                end
            else
                ind = s & flag==0;
                if sum(ind)>0
                    auto(tmp(ind)+1,sum(s),2) = auto(tmp(ind)+1,sum(s),2) + 1;
                end
            end
        end
    end
end

res = struct('autotime', autotime);
res = setfield(res, 'auto', auto);
head.NCounts = cnt;

