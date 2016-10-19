function [auto, rate, tcspcdata] = MultiFocus2FCS_Cluster_Calc(filename)
% [up down autotime] = MultiFocus2FCS_Cluster_Calc(filename)
% This is a helper for MultiFocus2FCS_Cluster that is invoked for each data
% package that is generated there

load(filename);

tcspcdata = zeros(NChannels,dnum);

if (~isempty(tmpy))
    for j=1:length(flv)
        dd = rev_dind(flv(j)+1);
        if (dd > 0)
            tcspcdata(tmpx(j)+1, dd) = tcspcdata(tmpx(j)+1, dd) + 1;
        end
    end
end    

flvp = zeros(length(flv),cnum*dnum);
r = 1;
for p=1:dnum
    for q=1:cnum
        flvp(:,r) = flv==dind(p) & tmpx>=tau(1,q,p) & max(tau(:,q,p))>=tmpx;
        r = r+1;
    end
end
rate = sum(flvp)/dt/Timeunit;

auto = tttr2xfcs(tmpy, flvp, Ncasc, Nsub);
end
