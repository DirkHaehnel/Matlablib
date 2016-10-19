function [err, c, z, c0, z0, zz] = AntibunchDelay0(p, t0, x, x0, bld)

% p(1) ... null time
% p(2) ... period
% p(3) ... delay
% p(4) ... decay time
% p(5) ... rotational diffusionc coefficients

t0 = t0(:)-p(1);
period = p(2);
t = abs(t0);
delay = mod(p(3),period);
rot = p(5);
tau = p(4);
n = 2*ceil(max(t)/period)*ceil(tau/period);

f = inline('sinh(x/tau).*exp(-delay/tau).*(x<=delay) + sinh(delay/tau).*exp(-x/tau).*(x>delay)','x','tau','delay');
ff = inline('exp((x<=delay).*(x-delay)/tau).*(x<=delay)/2 + exp(-(x>delay).*(x-delay)/tau).*(x>delay)/2','x','tau','delay');

if length(p)>4
    zz = [ones(length(t),1) x0(:)/mean(x0) zeros(length(t),6)];
    for j=0:n
        zz(:,3) = zz(:,3) + exp(-t/rot).*f(t,tau,j*period);
        zz(:,4) = zz(:,4) + exp(-10/3*t/rot).*f(t,tau,j*period);        
        zz(:,5) = zz(:,5) + exp(-t/rot).*f(t,tau,j*period+delay);
        zz(:,6) = zz(:,6) + exp(-10/3*t/rot).*f(t,tau,j*period+delay);
        zz(:,7) = zz(:,7) + exp(-t/rot).*f(t,tau,(j+1)*period-delay);
        zz(:,8) = zz(:,8) + exp(-10/3*t/rot).*f(t,tau,(j+1)*period-delay);
    end
else
    zz = [ones(length(t),1) zeros(length(t),3)];
end
% for j=-n:n
%     zz(:,2) = zz(:,2) + ff(t0,tau,j*period);
%     zz(:,3) = zz(:,3) + ff(t0,tau,j*period+delay);
%     zz(:,4) = zz(:,4) + ff(t0,tau,(j+1)*period-delay);
% end

if length(p)>4
%     c = [zz(:,1) x(:,2)/mean(x(:,2)) zz(:,5:7)]\x(:,1);
%     z = [zz(:,1) x(:,2)/mean(x(:,2)) zz(:,5:7)]*c;
    c = zz\x(:,1);
    z = zz*c;
    if nargin>3 && ~isempty(bld)
        plot(t0, x(:,1), 'o', t0, z);
        axis tight
        drawnow
    end
    err = sum(sum((x(:,1)-z).^2));
else
    c0 = lsqnonneg(zz,x);
    z0 = zz*c0;
    if nargin>3 && ~isempty(bld)
        plot(t0, x, 'o', t0, z0);
        axis tight
        drawnow
    end
    err = sum(sum((x-z0).^2)./t);
end

% close; para = Simplex('AntibunchDelay',[63 25 16 2 30 2],[0 24 0 0 0 0],[inf 26 inf inf inf 10],[],[],0.032*bin,y(:,10),weight,1);

