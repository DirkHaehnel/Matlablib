% Program for modeling double transits in SMD experiments
% flag 0 = reading the experimental data
% flag 1 = calculating the BSD for singles
% flag 2 = BHD by pure MCS
% flag 3 = calculating the TT-MF and TT-F distributions for singles
% flag 4 = calculating the TT-MF and TT-F distributions for doubles
% flag 5 = fitting the BSD to the experimental results

flag = 5
cd c:/joerg/doc/double

if flag == 0
	load distdata; dave = bsd(:,1:300)';
	t = 15:length(dave(:,1));
	pressure=[35 40 55 70 90 110 135 170 200];
	close all
	h = axes;
	plot(t, dave(t, 1));
	hold
	for k=2:9 plot(t, dave(t, k)/max(dave(t,k))*max(dave(t,1)) + (k-1)*0.001); end
	hold
	set(h, 'Ytick', []);
	for k =1:9
		text(260, 0.0004 + (k-1)*0.001, ['p = ', num2str(pressure(k))]);
	end
	xlabel('Burst size (photoelectrons)'); ylabel('Relative frequency'); times
	figure
	h = axes;
	t = 3:length(bhd(1,:));
	plot(t, bhd(1, t));
	hold
	for k=2:9 plot(t, bhd(k,t)/max(bhd(k,t))*max(bhd(1,t)) + (k-1)*0.0025); end
	hold
	set(h, 'Ytick', []);
	for k =1:9
		text(10.8, 0.0011 + (k-1)*0.0025, ['p = ', num2str(pressure(k))]);
	end
	xlabel('Burst height (photoelectrons)'); ylabel('Relative frequency'); times; drawnow
end

if flag==1 | flag==2
	close all
	load dave;
	dave = dave(:,1)';
	dave(71:length(dave)) = [];
	dave = dave/sum(dave);
	planck = 6.6262e-34; % Planck's constant in SI
	lightspeed = 2.9979e8; % speed of light in SI
	N = 50; % total number of time steps per simulation
	v = 2.4e4; % molecule's flow velocity in mum/sec
	phibl = 6e-5; % phi_bl
	sigma = 5e-7; % absorption cross section in mum^2
	I = 1.02e-3/planck/lightspeed*514e-9; % total laser power in photons/sec
	D = 43; % diffusion constant in mum^2/sec
	dt = 5e-5; % time step in sec
	ww = 10^2; % square of waist diameter of laser beam in mum
	% slit calculation
	% x-axis = axis of detection 
	% y-axis = axis of laser beam
	% z-axis = axis of capillary
	na = 0.85; % numerical aperture
	index = 1.33; % index of refraction
	psi = asin(na/index);
	k1 = 1/sin(psi);
	k2 = k1*cos(psi);
	d = 10; % slit width in object space in mum
	R = 5; % radius of delivery in mum
	L = 100; % distance between injection capillary and laser beam
	shiftx = 0; % x-shift in mum
	shifty = 0; % y-shift in mum
	z0 = -3*sqrt(ww);
	lam = 175; % hydrodynamic acceleration parameter
	told = 1; t0 = 0;
	while abs(t0-told)>1e-6
		told = t0;
		t0 = L/v + 1/lam*(1 - exp(-lam*t0))
	end
	t0 = t0 + z0/v;
	phif = 0.98*9e-3*sigma*2*I/pi/ww*dt; % Q_fl * eta * sigma * I * dt
	phi = 2*I/pi/ww*phibl*sigma*dt; % photobleaching factor
	tsh = 0.19; % Transit time threshold
	maxfluo = length(dave); % maximum number of detected photons
	M = 1:maxfluo;
	FAC = cumsum(log(M));
	L = 1:N;
	xi=sqrt(2*D*dt);
end

if flag==1
	bsm = zeros(1, maxfluo+1);
	for k=1:10000
		% random variable dr
		% with diffusion:
		x = xi*randn(1,N);
		y = xi*randn(1,N);
		z = xi*randn(1,N) + v*dt;
		% initial value
		r = R*sqrt(rand);
		theta = 2*pi*rand;
		x(1) = r*cos(theta) + sqrt(2*D*t0)*randn + shiftx;
		y(1) = r*sin(theta) + sqrt(2*D*t0)*randn + shifty;
		z(1) = z0 + sqrt(2*D*t0)*randn;
		% path calculation
		x = cumsum(x);
		y = cumsum(y);
		z = cumsum(z);
		% drawing the path:
		% r = sqrt(x.^2+y.^2) + i*z; plot(i*conj(r)); drawnow;
		% CEF
		tst1 = max([atan(-(y + d/2)./abs(x)); -psi*ones(size(y))]);
		tst2 = min([atan((d/2 - y)./abs(x)); psi*ones(size(y))]);
		delta1 = sqrt(1 - k1^2*sin(tst1).^2) + 1e-10;
		delta2 = sqrt(1 - k1^2*sin(tst2).^2) + 1e-10;
		res = k1*asin(k1*sin(tst2)) - k2*atan(k2*sin(tst2)./delta2) - k1*asin(k1*sin(tst1)) + k2*atan(k2*sin(tst1)./delta1);
		res = 2*sin(psi)/(2*pi*(1-cos(psi)))*res.*(imag(res)==0);
		% end CEF
		% bleaching and fluorescence
		fluo = phif*exp(-2/ww*(x.^2 + z.^2)).*res;
		cumfluo = round(cumsum(fluo));
		bb = phi*exp(-2*(x.^2+z.^2)/ww);
		bleach = [exp(-cumsum(bb(1:N-1))).*bb(2:N) exp(-sum(bb))];
		bleach = bleach/sum(bleach);
		% burst height and burst size
		for jj = 0:min(cumfluo(N), maxfluo)
			bsm(jj+1) = bsm(jj+1) + sum(bleach(cumfluo==jj));
		end
		if rem(k, 50)==0
			plot(bsm); drawnow;
		end
	end
%	calculation of burst size distribution
	Mi = 0.1:0.1:maxfluo;
	bsm =  interp1([0 M], bsm, Mi); bsm = bsm/sum(bsm);
	bs = zeros(1, maxfluo);
	for j=1:length(M) 
		bs(j) = bs(j) + exp(j*log(Mi) - Mi - FAC(j))*bsm;
	end
	bs = bs/sum(bs);
	save bsd bs
end

if flag==2
	load distdata; bhd = bhd(1,:);
	ff = 15; % max number of fluorescence photons per time bin
	bh = zeros(1,ff);
	for k=1:10000
		% random variable dr
		% with diffusion:
		x = xi*randn(1,N);
		y = xi*randn(1,N);
		z = xi*randn(1,N) + v*dt;
		% initial value
		r = R*sqrt(rand);
		theta = 2*pi*rand;
		x(1) = r*cos(theta) + sqrt(2*D*t0)*randn + shiftx;
		y(1) = r*sin(theta) + sqrt(2*D*t0)*randn + shifty;
		z(1) = z0 + sqrt(2*D*t0)*randn;
		% path calculation
		x = cumsum(x);
		y = cumsum(y);
		z = cumsum(z);
		% drawing the path:
		% r = sqrt(x.^2+y.^2) + i*z; plot(i*conj(r)); drawnow;
		% CEF
		tst1 = max([atan(-(y + d/2)./abs(x)); -psi*ones(size(y))]);
		tst2 = min([atan((d/2 - y)./abs(x)); psi*ones(size(y))]);
		delta1 = sqrt(1 - k1^2*sin(tst1).^2) + 1e-10;
		delta2 = sqrt(1 - k1^2*sin(tst2).^2) + 1e-10;
		res = k1*asin(k1*sin(tst2)) - k2*atan(k2*sin(tst2)./delta2) - k1*asin(k1*sin(tst1)) + k2*atan(k2*sin(tst1)./delta1);
		res = 2*sin(psi)/(2*pi*(1-cos(psi)))*res.*(imag(res)==0);
		% end CEF
		% bleaching and fluorescence, MCS
		fluo = phif*exp(-2/ww*(x.^2 + z.^2)).*res;
		fluo = cumsum(exp((0:ff)'*log(fluo+1e-10*(fluo==0)) - ones(ff+1,1)*fluo - [0 FAC(1:ff)]'*ones(1,N)));
		fluo = sum((ones(ff+1,1)*rand(1, N)) > fluo);
		bb = phi*exp(-2*(x.^2+z.^2)/ww);
		ind = sum(ones(1,N)*rand<exp(-cumsum(bb)));
		tmp = max(fluo(1:ind));
		if tmp>0
			bh(tmp) = bh(tmp) + 1;
		end
		if rem(k, 50)==0
			plot(1:ff,bh/max(bh),2:12,bhd(2:12)/max(bhd(2:12)))
			drawnow;
		end
	end
	save bhd bh
end

if flag==3 % BS-TT
	load trans1pe.his;
	trans = hist(trans1pe, 0:N);
	maxfluo1 = maxfluo + 1;
	scatter1 = zeros(N, maxfluo);
	ff = 12; % max number of fluorescence photons per time slice
	tt = 1:N;
	for k=1:100000
		% random variable dr
		% with diffusion:
		x = xi*randn(1,N);
		y = xi*randn(1,N);
		z = xi*randn(1,N) + v*dt;
		% initial value
		r = R*sqrt(rand);
		theta = 2*pi*rand;
		x(1) = r*cos(theta) + sqrt(2*D*t0)*randn + shiftx;
		y(1) = r*sin(theta) + sqrt(2*D*t0)*randn + shifty;
		z(1) = z0 + sqrt(2*D*t0)*randn;
		% path calculation
		x = cumsum(x);
		y = cumsum(y);
		z = cumsum(z);
		% CEF
		tst1 = max([atan(-(y + d/2)./absm(x)); -psi*ones(size(y))]);
		tst2 = min([atan((d/2 - y)./absm(x)); psi*ones(size(y))]);
		delta1 = sqrt(1 - k1^2*sin(tst1).^2) + 1e-10;
		delta2 = sqrt(1 - k1^2*sin(tst2).^2) + 1e-10;
		res = k1*asin(k1*sin(tst2)) - k2*atan(k2*sin(tst2)./delta2) - k1*asin(k1*sin(tst1)) + k2*atan(k2*sin(tst1)./delta1);
		res = 2*sin(psi)/(2*pi*(1-cos(psi)))*res.*(imag(res)==0);
		% end CEF
		% bleaching and fluorescence
		fluo = phif*exp(-2/ww*(x.^2 + z.^2)).*res;
		fluo = cumsum(exp((0:ff)'*log(fluo+1e-10*(fluo==0)) - ones(ff+1,1)*fluo - [0 FAC(1:ff)]'*ones(1,N)));
		bb = phi*exp(-2*(x.^2+z.^2)/ww);
		bleach = [exp(-cumsum(bb(1:N-1))).*bb(2:N) exp(-sum(bb))];
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
	end
	scatter1 = scatter1/sum(sum(scatter1));
	plot(M, dave, M, sum(scatter1)); title('dist1'); drawnow;
	figure
	h = pcolor(scatter1); set(h,'EdgeColor', 'none'); colormap(jet); title('scatter1');
	save ms1.mat scatter1
end

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

if flag==5 % bsd fitting
	load bsd;
	mx = 7;
	load distdata;
	dave = bsd(:,1:300)';
	len = length(dave(:,1));
	bs = [bs zeros(1,len - length(bs))];
	bb = bs;
	tmp = bs;
	for k = 2:mx
		tmp = conv(bs, tmp);
		tmp = tmp(1:len);
		bb = [bb; tmp];
	end
	t = 15:100;
	res = [];
	close all
	c = dave(t, 1)'/bs(t);
	h = axes;
	plot(t, c*bs(t), t, dave(t, 1), '.b');
	set(h, 'ytick', []);
	xlabel('Burst size (photoelectrons)'); ylabel('Relative frequency'); times; drawnow
	t = 15:len;
	for j = 1:9
		c = dave(t, j)'/bb(:,t);
		res = [res; c];
	end
	c = res;
	j = 9;
	figure
	h = axes;
	plot(1:len, c(j,:)*bb, t, dave(t, j), '.g');
	hold
	for k = 1:mx
		plot(1:len, c(j,k)*bb(k,:), 'b');
	end
	hold
	set(h, 'ytick', []);
	xlabel('Burst size (photoelectrons)'); ylabel('Relative frequency'); times; drawnow
	kappa = [];
	options = foptions;
	kappa = 0.5*ones(1,9);
	kappa = fmins('doubmin', kappa, options, [], c, mx);
	pressure=[35 40 55 70 90 110 135 170 200];
	figure
	fac = -log(kappa)/[ones(1,9); pressure];
	plot(pressure, -log(kappa),'o', pressure, fac*[ones(1,9); pressure], 'b');
	xlabel('Pressure (mm Hg)'); ylabel('Ratio tau_b/tau_m'); times
end

