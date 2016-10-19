function [err, c, z] = AntibunchRotoDiff(p, t0, x)

% p(1) ... null time
% p(2) ... period
% p(3) ... delay
% p(4) ... decay time
% p(5) ... rotational diffusionc coefficient

t0 = t0(:)-p(1);
period = p(2);
t = abs(t0);
delay = mod(p(3),period);
tau = p(4);
rot = 1/p(5);
n = 2*ceil(max(t)/period)*ceil(tau/period);

f = inline('sinh(x/tau).*exp(-delay/tau).*(x<=delay) + sinh(delay/tau).*exp(-x/tau).*(x>delay)','x','tau','delay');
ff = inline('exp((x<=delay).*(x-delay)/tau).*(x<=delay)/2 + exp(-(x>delay).*(x-delay)/tau).*(x>delay)/2','x','tau','delay');

tmp = exp(-6*t*rot(1));
zz = zeros(length(t),2); % ones(length(t),1)
for j=0:n
    zz(:,2) = zz(:,2) + tmp.*f(t,tau,j*period+delay);
end
for j=-n:n
    zz(:,1) = zz(:,1) + ff(t0,tau,j*period+delay);
end

c = zz\x;
z = zz*c;
plot(t0, x, 'o', t0, z);
axis tight
drawnow
err = sum(sum((x-z).^2));

% close; para = Simplex('AntibunchDelay',[63 25 16 2 30 2],[0 24 0 0 0 0],[inf 26 inf inf inf 10],[],[],0.032*bin,y(:,10),weight,1);

