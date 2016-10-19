function [err, c1, c2, z, zz] = AntibunchDelayOld(p, t0, x)

% p(1) ... null time
% p(2) ... period
% p(3) ... delay
% p(4) ... decay time
% p(5) ... rotation time
% p(6) ... bleed-through in detection

t0 = t0(:)-p(1);
period = p(2);
t = abs(t0);
delay = mod(p(3),period);
tau = p(4);
rot = p(5);
if nargin<5 || isempty(las) 
    las = 1;
end
n = 2*ceil(max(t)/period)*ceil(tau/period);

f = inline('sinh(x/tau).*exp(-delay/tau).*(x<=delay) + sinh(delay/tau).*exp(-x/tau).*(x>delay)','x','tau','delay');
ff = inline('exp((x-delay)/tau).*(x<=delay)/2 + exp(-(x-delay)/tau).*(x>delay)/2','x','tau','delay');

weight = zeros(3,size(aniso.ux,1));
weight(1,1) = real(aniso.ux(1,1).*conj(aniso.uy(1,1)))/2;
weight(2,1) = real(aniso.ux(1,1).*conj(aniso.uyr(1,1)))/2;
weight(3,1) = real(aniso.uy(1,1).*conj(aniso.uxr(1,1)))/2;
for L = 2:size(aniso.ux,1)
    for M = 1:L
        weight(1,L) = weight(1,L)+ real(aniso.ux(L,M).*conj(aniso.uy(L,M)));
        weight(2,L) = weight(1,L)+ real(aniso.ux(L,M).*conj(aniso.uyr(L,M)));
        weight(3,L) = weight(1,L)+ real(aniso.uy(L,M).*conj(aniso.uxr(L,M)));
    end
end

weight = diag([1+las^2 las las 1+las^2 las])*weight;

for L = 1:size(aniso.ux,1)
    epp = epp + weight(1,L)*exp(-L*(L-1)*t/6/rot);
    epo = epo + weight(2,L)*exp(-L*(L-1)*t/6/rot);
    eop = eop + weight(3,L)*exp(-L*(L-1)*t/6/rot);
end

zz = zeros(length(t),2);
for j=-n:n
    zz(:,1) = zz(:,1) + weight(1,1)*ff(t0,tau,j*period) + ...
        weight(2,1)*ff(t0,tau,j*period+delay) + ...
        weight(3,1)*ff(t0,tau,(j+1)*period-delay);
end
for j=0:n
    zz(:,2) = zz(:,2) + epp.*f(t,tau,j*period) + ...
        (epo.*(t0>=0) + eop.*(t0<0)).*f(t,tau,j*period+delay) + ...
        (eop.*(t0>=0) + epo.*(t0<0)).*f(t,tau,(j+1)*period-delay);
end

zz = [ones(size(x,1),1) zz];
c1 = lsqnonneg(zz(t0>=0,:),x(t0>=0));
c2 = lsqnonneg(zz(t0<0,:),x(t0<0));
z = [zz(t0<0,:)*c2; zz(t0>=0,:)*c1];

if nargin>5 && ~isempty(bld)
    plot(t0, x, 'o', t0, z); 
    drawnow
end

err = sum(sum((x-z).^2)./z./t);
%err = sum(sum((x-z).^2))

% close; para = Simplex('AntibunchDelay',[63 25 16 2 30 2],[0 24 0 0 0 0],[inf 26 inf inf inf 10],[],[],0.032*bin,y(:,10),weight,1);

