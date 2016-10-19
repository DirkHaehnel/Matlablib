% Programm zur Simulation eines FCS-Experiments

close all

a = 0.5;
b = 2;
N = 2^17; % max number of time channels
dt = 1e-6; % time interval in microseconds
L = 1:N;
time = (1:N)*dt;
dif = 1e3; % diffusion constant in mum^2/time unit
auto = zeros(1,N);
xi = sqrt(2*dif*dt);
rhomax = 2.5*a;
zmax = 2.5*b;
fluototal = 0;
x = 0;
y = 0;
z = 0;
for k=1:1000
   % random path
   x = x(1) + cumsum(xi*randn(1,N));
   y = y(1) + cumsum(xi*randn(1,N));    
   z = z(1) + cumsum(xi*randn(1,N));
    
   % fluorescence
   rr(rr>rhomax) = rhomax;
   zz(zz>zmax) = zmax;
   rr(rr<rhomin) = rhomin;
   zz(zz<zmin) = zmin;
   fluo = interp2(z, rho, f0, zz, rr) + dt*background;
   fluo = cumsum(cumprod(repmat(fluo,[3 1])));
   fluo = sum(ones(3,1)*(rand(1,N).*fluo(end,:))>fluo);
   fluototal = fluototal + sum(fluo);
   flc = fft([fluo zeros(1,N)]);
   flc = abs(ifft(flc.*conj(flc)));
   auto = auto + flc(1:N);
   tmp = auto/k/N - (fluototal/N/k)^2;
   semilogx(time(2:end), auto(2:end)); drawnow;
end

