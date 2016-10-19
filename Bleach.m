% Programm zur Simulation eines FCS-Experiments

N = 2^11; % Gesamtzahl der Zeitkanäle
D = 43; % Diffusionskonstante in mum^2/sec
dt = 3e-3; % Zeitintervall in sec
ar = 25; % halbe Kastenkante a
hr = 25; % halbe Kastenkante h
L = 1:N;
a2D = 1.6791;
c2D = 0.1452;
w2D = 2;
a2 = a2D*2*D;
c2 = c2D*2*D;
w2 = w2D*2*D;
tt=(1:N)*dt; 
tmpa = a2D + tt;
tmpc = c2D + tt;
a0 = a2D*sqrt(c2D)./tmpa./sqrt(tmpc);
phi = 57; % bleaching parameter
auto = phi*sqrt(tmpa./(tmpa+w2D)).*atanh(tt./sqrt(tmpa.*(tmpa+w2D)));
auto = a0.*exp(-auto);
phi = phi/w2D*dt;
close all
semilogx(tt,a0/a0(1),tt,auto/auto(1)); drawnow; hold on;
bild = semilogx(tt,zeros(size(tt)), 'EraseMode', 'xor');

xi=sqrt(2*D*dt);
auto = zeros(1,N);
mm = 0;

for k=1:10000
	% random variable dr
	x = xi*randn(1,N);
	y = xi*randn(1,N);
	z = xi*randn(1,N);
	% initial value
	x(1) = -ar; % ar*(2*rand(1)-1);
	y(1) = ar*(2*rand(1)-1);
	z(1) = hr*(2*rand(1)-1);
	% path calculation
	x = cumsum(x);
	y = cumsum(y);
	z = cumsum(z);
	r = sqrt(x.^2+y.^2) + i*z;
%	plot(r); drawnow;
	% fluorescence
	fluo = exp(-real(r).^2/a2-imag(r).^2/c2);
	t = L(rand>exp(-phi*cumsum(exp(-2*real(r).^2/w2))));
	fluo(t) = zeros(size(t));
	mm = mm + sum(fluo);
%	plot(fluo); drawnow;
	flc = fft([fluo zeros(1,N)]);
	flc = real(ifft(flc.*conj(flc)));
	auto = auto + flc(1:N);
	tmp = auto/k/N - (mm/N/k)^2;
	if rem(k,10)==0 set(bild, 'ydata', tmp/tmp(1)); drawnow; end;
end

% save
auto = auto/k/N;
mm = mm/k/N;
% save bleach.mat auto mm f1 f2 bl D dt N phi ar hr

% laser intensity
% int=exp(-2*rho^2/(w0^2+(z*lam/pi/w0^2)^2))*(2/(w0^2+(z*lam/pi/w0^2)^2)/pi);
