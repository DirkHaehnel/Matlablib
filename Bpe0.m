% Programm for the simulation of Phycoerythrin crossings through a laser beam
% adapted for simulating the inlfuence of gackground
% no diffusion, no photobleaching

close
N = 100; % total number of time steps per simulation
v = 1e3; % molecule's flow velocity in mum/sec
dt = 2.5e-5; % time step in sec
ww = 10^2; % square of waist diameter of laser beam in mum
% slit calculation
% x-axis = axis of detection 
% y-axis = axis of laser beam
% z-axis = axis of capillary
L = 100; % distance between injection capillary and laser beam
z0 = -3*sqrt(ww);
told = 1; t0 = 0;
t0 = z0/v;
phif = 5e6;
bg = 10/N; % background
maxfluo = 100; % maximum number of detected photons
M = 1:maxfluo;
FAC = cumsum(log(M));
L = 1:N;

xi=sqrt(2*D*dt);
mm = [];
dist = zeros(1,maxfluo);

for k=1:1
	% random variable dr
	% with diffusion:
	x = zeros(1,N);
	y = zeros(1,N);
	z = v*dt*ones(1,N);
	% initial value
	z(1) = z0;
	% path calculation
	z = cumsum(z);
	% fluorescence
	fluo = bg + phif*exp(-2/ww*(x.^2 + z.^2));
	% bleaching sampling
	mm = sum(fluo);
	dist = dist + exp(log(mm)*M - mm - FAC);
	plot(M, dist/sum(dist)); title(['path # = ', num2str(k)]);
end

% exact laser intensity
% int=exp(-2*rho^2/(w0^2+(z*lam/pi/w0^2)^2))*(2/(w0^2+(z*lam/pi/w0^2)^2)/pi);
