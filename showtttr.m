function showtttr(fname, binwidth);

if nargin<2
    binwidth = 1e3;
end

photons = 1e4; % number of photons processed one at a time
close all   

[tmp, tmp, tmp, head] = tttrRead(fname,[1 1]);
timeunit = head.GlobClock*1e-9;    

tst = 1;
time = 0;
cnt = 0;

ax1 = [0 0 0 0]; ax2 = [0 0 0 0];
while tst>0
    [x, flag, tmp, head, tst, overcount, marker] = tttrRead(fname,[cnt+1 photons]);        
    x(marker>0) = [];
    flag(marker>0) = [];
    if tst>0
        cnt = cnt+tst;
        y = tttr2bin(x,binwidth,flag);
        subplot(211);
        t = (time+(1:length(y))*binwidth)*1e-7;
        plot(t,y(:,1));
        tmp = axis;
        if tmp(4)>ax1(4)
            ax1(4) = tmp(4);
        end
        axis([t([1 end]) ax1(3:4)]);
        subplot(212);
        plot(t,y(:,2));        
        tmp = axis;
        if tmp(4)>ax2(4)
            ax2(4) = tmp(4);
        end
        axis([t([1 end]) ax2(3:4)]);        
        drawnow
        dt = x(end) - x(1);
        time = time + dt;
    end
end

