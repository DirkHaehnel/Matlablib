function res = PhaseMod2FCS(fname, period, maxtime)

photons = 1e5; % number of photons processed one at a time

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

if strcmp(fname(end-2:end),'spc')
    tst=1;
    jj=1;
    [head, head, head, head] = Spcread_4096(fname);
    Timeunit = head.GlobClock*1e-9;
    while tst
        [y, flag, num, num, num] = Spcread_4096(fname, [cnt+1 photons]);
%         auto1(:,jj) = tttr2xfcs(floor(y*Timeunit/period),[flag==0 flag==1],Ncasc,Nsub);
%         auto2(:,jj) = abs(tttr2xfcs(floor(y*Timeunit/period),[flag==0 flag==1].*(exp(2*pi*i/period*y*Timeunit)*[1 1]),Ncasc,Nsub));
        auto1 = auto1 + tttr2xfcs(floor(y*Timeunit/period),[flag==0 flag==1],Ncasc,Nsub);
        auto2 = auto2 + abs(tttr2xfcs(floor(y*Timeunit/period),[flag==0 flag==1].*(exp(2*pi*i/period*y*Timeunit)*[1 1]),Ncasc,Nsub));
        time(jj) = dt*Timeunit;
        semilogx(autotime, sum(auto1+4*auto2,2)/sum(time),autotime, sum(auto1-4*auto2,2)/sum(time));
        drawnow;
        jj = jj+1;
        if num<photons
            tst = 0;
        end
    end
else
    Timeunit = 12.5e-9;
    fin = fopen(fname);
    fgetl(fin);
    while tst>0
        x = fscanf(fin,'%d',2*photons);
        tst = length(x);
        if tst>0
            x = reshape(x,2,length(x)/2)';
            ind = x<[0 0;x(1:end-1,:)];
            x = x + cumsum(ind*(2^32-1));
            cnt = cnt + tst;
            dt = x(end,1)-x(1,1);
            if dt>2*autotime(end)
                time(jj) = dt*Timeunit;

                auto1(:,jj) = tttr2xfcs(x(:,1)*Timeunit/period, ones(length(x(:,1)),1), Ncasc, Nsub);
                auto2(:,jj) = abs(tttr2xfcs(round(x(:,1)*Timeunit/period), exp(2*pi*i/period*x(:,1)*Timeunit), Ncasc, Nsub));

                semilogx(autotime, sum(auto1+4*auto2,2)/sum(time),autotime, sum(auto1-4*auto2,2)/sum(time));
                semilogx(autotime, sum(auto1,2)/sum(time));
                drawnow;
            end
        end
        jj = jj+1;
    end
    fclose(fin);
end

for j=1:length(time)
    auto1(:,j) = auto1(:,j)/time(j)/Timeunit;
%    auto2(:,j) = auto2(:,j)/time(j)/Timeunit;
end

clf
semilogx(autotime, sum(auto1+4*auto2,2),autotime, sum(auto1-4*auto2,2));
res.autotime = autotime;
res.auto1 = auto1;
% res.auto2 = auto2;
% res.auto = auto1+4*auto2;
% res.cross = auto1-4*auto2;
res.rate = rate;
res.time = time;
res.NCounts = cnt;


