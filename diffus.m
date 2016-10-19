% Programm zur Simulation eines FCS-Experiments

N = 2^12; % Gesamtzahl der Zeitkanäle
D = 10^3; % Diffusionskonstante in mum^2/sec
dt = 10^(-5); % Zeitintervall in sec
f1 = 4; % Beobachtungsfleck - grosse Halbachse in mum
f2 = 2; % Beobachtungsfleck - kleine Halbachse in mum
bl = 4; % Bleachingfleck - grosse Halbachse in mum
phi = 0; % phi_bl/phi_f
ar = 5; % halbe Kastenkante a
hr = 5; % halbe Kastenkante h
L = 1:N;
tt=(1:N)*dt; tst=1./(f1^2+4*D*tt)./sqrt(f2^2+4*D*tt); tst=tst/tst(1);

xell = 0:0.1:f2; yell = f1*sqrt(1-xell.*xell/f2^2); xell = [xell fliplr(xell)]; yell = [yell -fliplr(yell)];
semilogx(tt,tst); drawnow; hold on;
bild = semilogx(tt,zeros(size(tt)), 'EraseMode', 'xor');

xi=sqrt(2*D*dt);
auto = zeros(1,N);
mm = 0;
for k=1:1000
	% random variable dr
	x = xi*randn(1,N);
	y = xi*randn(1,N);
	z = xi*randn(1,N);
	% initial value
	x(1) = ar*(2*rand(1)-1);
	y(1) = ar*(2*rand(1)-1);
	z(1) = hr*(2*rand(1)-1);
	% path calculation
	x = cumsum(x);
	y = cumsum(y);
	z = cumsum(z);
	r = sqrt(x.^2+y.^2) + i*z;
%	plot(xell, yell); axis([0 30 -20 20]); hold;
%	plot(r); hold; drawnow;
	% fluorescence
	fluo = exp(-2*real(r).^2/f1^2-2*imag(r).^2/f2^2);
	mm = mm + sum(fluo);
%	t = L(rand<exp(-phi*cumsum(exp(-2*real(r).^2/bl^2))));
%	fluo(t) = zeros(size(t));
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
% save diffus.mat auto mm f1 f2 bl D dt N phi ar hr

% laser intensity
% int=exp(-2*rho^2/(w0^2+(z*lam/pi/w0^2)^2))*(2/(w0^2+(z*lam/pi/w0^2)^2)/pi);
% tt=(1:N)*dt; tst=1./(f1^2+4*D*tt)./sqrt(f2^2+4*D*tt); tst=tst/tst(1);
% semilogx(tt,tst,tt,auto/auto(1))
% cr = phi*bl^2/2/D/dt*sqrt((f1^2+8*D*tt)./(f1^2+4*bl^2+8*D*tt)).*atanh(8*D*tt./sqrt((f1^2+8*D*tt).*(f1^2+4*bl^2+8*D*tt)));

