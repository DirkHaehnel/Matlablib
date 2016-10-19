function res = BurstSelect(res,sel)

tt = 1:length(res.sync);
t1 = res.t1(sel);
t2 = res.t2(sel);
ind1 = zeros(numel(res.sync)+1,1);
ind2 = ind1;
ind1(tt(t1)) = 1;
ind2(tt(t2)+1) = 1;
ind1 = (cumsum(ind1)-cumsum(ind2)>0).*cumsum(ind1);
ind1(end)=[];
res.sync = res.sync(ind1>0);
res.tcspc = res.tcspc(ind1>0);
res.chan = res.chan(ind1>0);
res.ind = res.ind(ind1>0);
res.bs = res.bs(sel,:);
ind1 = ind1(ind1>0);
tt = 1:numel(ind1);
res.t1 = tt(diff([0; ind1(1:end)])==1);
res.t2 = tt(diff([ind1; max(ind1)+1])==1);
    

