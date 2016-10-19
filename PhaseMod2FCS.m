function [res, head] = PhaseMod2FCS(fname, period, maxtime)

photons = 5e5; % number of photons processed one at a time

close all

if nargin<2 || isempty(period)
    period = 10e-6;
end

if nargin<3 || isempty(maxtime)
    Nsub = 10;
    Ncasc = ceil(log2(3/period/Nsub)); % 3 sec max correlation time
else
    if length(maxtime)==2
        Nsub = maxtime(2);
    else
        Nsub = 10;
    end
    Ncasc = ceil(log2(maxtime(1)/period/Nsub));
end
autotime = cumsum(reshape(repmat(2.^(0:Ncasc-1),Nsub,1),Ncasc*Nsub,1));
autotime = autotime*period;
auto1 = zeros(length(autotime),2,2);
auto2 = auto1;

rate = [];
time = 0;
cnt = 0;
tst = 1;
jj = 1;
[head, head, head, head] = Spcread_4096(fname,[0 0]);
Timeunit = head.GlobClock*1e-9;
while tst
    [y, flag, num, num, num] = Spcread_4096(fname, [cnt+1 photons]);
    dt = y(end)-y(1);
    cnt = cnt + num;
    auto1(:,:,:,jj) = dt*tttr2xfcs(floor(y*Timeunit/period),[flag==0 flag==1],Ncasc,Nsub);
    auto2(:,:,:,jj) = dt*abs(tttr2xfcs(floor(y*Timeunit/period),[flag==0 flag==1].*(exp(2*pi*i/period*y*Timeunit)*[1 1]),Ncasc,Nsub));
%     auto1(:,:,:,jj) = dt*tttr2xfcs(y,[flag==0 flag==1],Ncasc,Nsub);
%     auto2(:,:,:,jj) = dt*abs(tttr2xfcs(y,[flag==0 flag==1].*(exp(2*pi*i/period*y*Timeunit)*[1 1]),Ncasc,Nsub));
    time(jj) = dt*Timeunit;
    rate(jj) = num;
    semilogx(autotime, squeeze(sum(auto1(:,1,2,:)+auto1(:,2,1,:)+4*auto2(:,1,2,:)+4*auto2(:,2,1,:),4)),...
        autotime, squeeze(sum(auto1(:,1,2,:)+auto1(:,2,1,:)-4*auto2(:,1,2,:)-4*auto2(:,2,1,:),4)));
    drawnow;
    jj = jj+1;
    if num<photons
        tst = 0;
    end
end

auto1 = auto1*Timeunit/sum(time);
auto2 = auto2*Timeunit/sum(time);

res.autotime = autotime;
res.auto1 = auto1;
res.auto2 = auto2;
res.auto = auto1+4*auto2;
res.cross = auto1-4*auto2;
res.rate = rate;
res.time = time;
res.NCounts = cnt;


