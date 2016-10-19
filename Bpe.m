% Programm for the simulation of Phycoerythrin crossings through a laser beam

close
planck = 6.6262e-34; % Planck's constant in SI
lightspeed = 2.9979e8; % speed of light in SI
N = 100; % total number of time steps per simulation
v = 2.4e4; % molecule's flow velocity in mum/sec
phibl = 6e-5; % photodestruction rate
sigma = 5e-7; % absorption cross section in mum^2
lambda = 514e-9; % excitation wavelength in m
I = 1.02e-3/planck/lightspeed*lambda; % total laser power in photons/sec
D = 43; % diffusion constant in mum^2/sec
dt = 2.5e-5; % time step in sec
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
R = 0; % radius of delivery in mum
L = 100; % distance between injection capillary and laser beam
shiftx = 0; % x-shift in mum = x-displacement of molecule flow from z-axis
shifty = 0; % y-shift in mum = y-displacement of molecule flow from z-axis
z0 = -3*sqrt(ww); % start plane of molecules
lam = 175; % empirical hydrodynamic acceleration parameter
told = 1; t0 = 0;
while abs(t0-told)>1e-6
	told = t0;
	t0 = L/v + 1/lam*(1 - exp(-lam*t0))
end
t0 = t0 + z0/v;
phif = 0.98*9e-3*sigma*2*I/pi/ww*dt; % fluorescence quantum yield * overall detection efficiency * sigma * I * dt
phi = 2*I/pi/ww*phibl*sigma*dt; % photobleaching factor
bg = 100; % total background within the time window
maxfluo = 100; % maximum number of detected photons
M = 1:maxfluo;
FAC = cumsum(log(M));
L = 1:N;

xi=sqrt(2*D*dt);
mm = [];
dist = zeros(1,maxfluo);

for k=1:1000
	% random variable dr
	% with diffusion:
	x = xi*randn(1,N);
	y = xi*randn(1,N);
	z = xi*randn(1,N) + v*dt;
	% initial position
	r = R*sqrt(rand);
	theta = 2*pi*rand;
	x(1) = r*cos(theta) + sqrt(2*D*t0)*randn + shiftx;
	y(1) = r*sin(theta) + sqrt(2*D*t0)*randn + shifty;
	z(1) = z0;
	% path calculation
	x = cumsum(x);
	y = cumsum(y);
	z = cumsum(z);
	% drawing the path:
%	r = sqrt(x.^2+y.^2) + i*z; plot(i*conj(r)); drawnow;
	% fluorescence
	% CEF
	tst1 = max([atan(-(y + d/2)./abs(x)); -psi*ones(size(y))]);
	tst2 = min([atan((d/2 - y)./abs(x)); psi*ones(size(y))]);
	delta1 = sqrt(1 - k1^2*sin(tst1).^2) + 1e-10;
	delta2 = sqrt(1 - k1^2*sin(tst2).^2) + 1e-10;
	res = k1*asin(k1*sin(tst2)) - k2*atan(k2*sin(tst2)./delta2) - k1*asin(k1*sin(tst1)) + k2*atan(k2*sin(tst1)./delta1);
	res = 2*sin(psi)/(2*pi*(1-cos(psi)))*res.*(imag(res)==0);
	% end CEF
	fluo = bg + phif*exp(-2/ww*(x.^2 + z.^2)).*res;
	% bleaching sampling
	mm = sum(fluo);
	if mm>0
		bb = phi*exp(-2*(x.^2+z.^2)/ww);
		bleach = exp(-cumsum(bb(1:N-1))).*bb(2:N);
		mm = cumsum(fluo);
		for kk = 1:N-1
			if mm(kk)>0 
				dist=dist+bleach(kk)*exp(log(mm(kk))*M-mm(kk)-FAC); 
			end
		end
		dist = dist + exp(-sum(bb))*exp(log(mm(N))*M -mm(N)-FAC);
	end
	if rem(k, 50)==0
		plot(M, dave/sum(dave), M, dist/sum(dist)); title(['path # = ', num2str(k)]);
		drawnow;
	end
end

% exact laser intensity distribution
% int=exp(-2*rho^2/(w0^2+(z*lam/pi/w0^2)^2))*(2/(w0^2+(z*lam/pi/w0^2)^2)/pi);
