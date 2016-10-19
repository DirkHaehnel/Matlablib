% Program for modeling transits in 2f-FCS-setup

close all
clear all
NatConst

Nmax = 1e6;
lamex = 0.64;
lamem = 0.67;
delta = 0.2/sqrt(2);
wLaser = 0.4;
wDetect = 0.5;
av = 100/60;
v = 1e4; % molecule's flow velocity in mum/sec
sigma = 1e-8; % absorption cross section in mum^2
power = 100e-6/PlanckConstant/SpeedOfLight*lamex*1e-6; % total laser power in photons/sec
D = 1e1; % diffusion constant in mum^2/sec
dt = 5e-6; % time step in sec
n = 1.33; % index of refraction
z0 = -10; % start point
x0 = 10;
N = round(2*abs(z0)/v/dt); % total number of time steps per simulation
L = 1:N;
eta = 0.1;
qy = 0.5;
phif = qy*eta*sigma*2*power/pi*dt; % Q_fl * eta * sigma * I * dt

xi=sqrt(2*D*dt);

tdist = zeros(1,Nmax); t1 = tdist; t2 = tdist; bs1 = tdist; bs2 = tdist; cnt = 1;
for k=1:Nmax
    % path calculation
    x = x0*(2*rand-1) + cumsum(xi*randn(1,N));
    y = x0*(2*rand-1) + cumsum(xi*randn(1,N));
    z = z0 + cumsum(xi*randn(1,N) + v*dt);
    % CEF
    res = D1CEF([wDetect 0],av,lamem,x);
    % fluorescence
    ww = LaserBeamFit([wLaser 0],lamex,x).^2;
    fluo1 = Poissrnd(phif*res.*exp(-2*((y+delta).^2 + (z+delta).^2)./ww)./ww);
    fluo2 = Poissrnd(phif*res.*exp(-2*((y-delta).^2 + (z-delta).^2)./ww)./ww);
    if sum(fluo1)>=10 && sum(fluo2)>=10
        plot(1:N,fluo1,1:N,fluo2); drawnow
        tdist(cnt) = (sum(L.*fluo2)./sum(fluo2)-sum(L.*fluo1)./sum(fluo1)).*dt;
        t1(cnt) = sqrt(sum(L.^2.*fluo1)./sum(fluo1)-(sum(L.*fluo1)./sum(fluo1))^2).*dt;
        t2(cnt) = sqrt(sum(L.^2.*fluo2)./sum(fluo2)-(sum(L.*fluo2)./sum(fluo2))^2).*dt;
        if rand>0.5
            bs1(cnt) = 2*sum(fluo1);
            bs2(cnt) = 2*sum(fluo2);
        else
            bs1(cnt) = sum(fluo1);
            bs2(cnt) = sum(fluo2);
        end
        cnt = cnt+1;
    end
    if rem(k, Nmax/1e2)==0
        mhist(tdist(1:cnt-1)); drawnow;
    end
end
tdist(cnt:end) = []; t1(cnt:end) = []; t2(cnt:end) = []; bs1(cnt:end) = []; bs2(cnt:end) = []; 

plot3(log(abs(t1./tdist)),log(bs1./bs2),bs1+bs2,'o')

return

ind=t1<=1.5*tdist & t2<=1.5*tdist & 0.5*tdist<=t1 & 0.5*tdist<=t2; 
a = log(bs1(ind)./bs2(ind));
b = bs1(ind)+bs2(ind);
c = abs(t1(ind)./tdist(ind));
[a,ord]=sort(a);
b=b(ord);
c=c(ord);
close; pg = simplex('Gauss',[0 1],[0 0],[0 inf],[],[],a,b,[],1,1);

weight = exp(-(a-pg(1)).^2/pg(2)^2/2);
mhist(b./weight,(0.5:100)/100*3*max(b),weight.^5)
[y,x] = mhist(b./weight,(0.5:100)/100*1.2*max(b),weight.^5);
close; pg = simplex('Gauss',[100 200 10 10],[],[],[],[],x,y,[],1,1);


for kk = 2:N-1
    ind = sum((ones(ff+1,1)*rand(1,kk)) > fluo(:,1:kk));
    t2 = tt([([ind(2:kk) 0]==0 & ind>0) zeros(1,N-kk)]);
    t1 = tt([([0 ind(1:kk-1)]==0 & ind>0) zeros(1, N-kk)]);
    transit = t2-t1;
    for jj=1:length(t1)
        ind = sum(ind(t1(jj):t2(jj)));
        if (ind>0 & ind<maxfluo1 & transit(jj)>0)
            scatter1(transit(jj), ind) = scatter1(transit(jj), ind) + bleach(kk);
        end
    end
end
ind = sum((ones(ff+1,1)*rand(1,N)) > fluo);
t2 = tt([ind(2:N) 0]==0 & ind>0);
t1 = tt([0 ind(1:N-1)]==0 & ind>0);
transit = t2-t1;
for jj = 1:length(t1)
    ind = sum(ind(t1(jj):t2(jj)));
    if (ind>0 & ind<maxfluo1 & transit(jj)>0)
        scatter1(transit(jj), ind) = scatter1(transit(jj), ind) + bleach(N);
    end
end
if rem(k, 50)==0
    %			plot(0:N, trans/sum(trans), 1:N, sum(scatter1')/sum(sum(scatter1')))
    h = pcolor(scatter1(1:size(scatter1)/2,1:size(scatter1')/1.5)); set(h,'EdgeColor', 'none'); colormap(jet);
    drawnow;
end
k

if flag==4 % BS-TT for multiple events
	load ms1;
	maxfluo = 120;
	M = 1:maxfluo;
	scatter2 = zeros(N+1, maxfluo+1);
	ind = 1:N+1;
	ind = ind(sum(scatter1')>0);
	scatterf = fft([scatter1 zeros(N+1, maxfluo - size(scatter1'))]')';
	for j = 1:length(ind)
		for k = j:length(ind)
			tmp = scatterf(ind(j),:).*scatterf(ind(k),:);
			scatter2(ind(k),:) = scatter2(ind(k),:) + (ind(k)-ind(j)+1)*tmp;
			for kk = 1:ind(j)
				scatter2(ind(k)+kk,:) = scatter2(ind(k)+kk,:) + 2*tmp;
			end
		end
	end
	scatter2 = absm(ifft(scatter2'))';
	scatter2 = scatter2/sum(sum(scatter2));
	figure
	plot(1:length(sum(scatter1)), sum(scatter1), M, sum(scatter2)); drawnow;
	figure
	h = pcolor(scatter2); set(h,'EdgeColor', 'none'); colormap(jet); title('scatter2');
	save ms2.mat scatter2
end
