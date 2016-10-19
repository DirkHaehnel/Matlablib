function [err, c, z, c0, z0, zz] = AntibunchDelay(p, t0, x, aniso, las, bld)

% p(1) ... null time
% p(2) ... period
% p(3) ... delay
% p(4) ... decay time
% p(5) ... depolarization angle
% p(6:7) ... rotational diffusionc coefficients

t0 = t0(:)-p(1);
period = p(2);
t = abs(t0);
delay = mod(p(3),period);
tau = p(4);
psi = p(5);
if length(p)>5
    if length(p)==6
        rot = [p(6) 0];
    else
        rot = p(6:7);
    end
end
if nargin<5 || isempty(las) 
    las = 1;
end
n = 2*ceil(max(t)/period)*ceil(tau/period);

f = inline('sinh(x/tau).*exp(-delay/tau).*(x<=delay) + sinh(delay/tau).*exp(-x/tau).*(x>delay)','x','tau','delay');
ff = inline('exp((x<=delay).*(x-delay)/tau).*(x<=delay)/2 + exp(-(x>delay).*(x-delay)/tau).*(x>delay)/2','x','tau','delay');

c = cos(psi);
s = sin(psi);

weight = zeros(3,size(aniso.ux,1),size(aniso.ux,2));
weight(1,1,1) = real((c^4*aniso.ux(1,1)+c^2*s^2*(aniso.uxr(1,1)+aniso.uy(1,1))+s^4*aniso.uyr(1,1))*...
        (c^4*aniso.uy(1,1)+c^2*s^2*(aniso.ux(1,1)+aniso.uyr(1,1))+s^4*aniso.uxr(1,1)) + ...
        las^2*(c^4*aniso.uxr(1,1)+c^2*s^2*(aniso.uyr(1,1)+aniso.ux(1,1))+s^4*aniso.uy(1,1))*...
        (c^4*aniso.uyr(1,1)+c^2*s^2*(aniso.uxr(1,1)+aniso.uy(1,1))+s^4*aniso.ux(1,1)))/2;
weight(2,1,1) = las*real((c^4*aniso.ux(1,1)+c^2*s^2*(aniso.uxr(1,1)+aniso.uy(1,1))+s^4*aniso.uyr(1,1))*...
        conj(c^4*aniso.uyr(1,1)+c^2*s^2*(aniso.uxr(1,1)+aniso.uy(1,1))+s^4*aniso.ux(1,1)))/2;
weight(3,1,1) = las*real((c^4*aniso.uy(1,1)+c^2*s^2*(aniso.uyr(1,1)+aniso.ux(1,1))+s^4*aniso.uxr(1,1))*...
        (c^4*aniso.uxr(1,1)+c^2*s^2*(aniso.ux(1,1)+aniso.uyr(1,1))+s^4*aniso.uy(1,1)))/2;
for L = 2:size(aniso.ux,1)
    M = 1;
    weight(1,L,M) = real((c^4*aniso.ux(L,M)+c^2*s^2*(aniso.uxr(L,M)+aniso.uy(L,M))+s^4*aniso.uyr(L,M))*...
        (c^4*aniso.uy(L,M)+c^2*s^2*(aniso.ux(L,M)+aniso.uyr(L,M))+s^4*aniso.uxr(L,M)) + ...
        las^2*(c^4*aniso.uxr(L,M)+c^2*s^2*(aniso.uyr(L,M)+aniso.ux(L,M))+s^4*aniso.uy(L,M))*...
        (c^4*aniso.uyr(L,M)+c^2*s^2*(aniso.uxr(L,M)+aniso.uy(L,M))+s^4*aniso.ux(L,M)))/2;
    weight(2,L,M) = las*real((c^4*aniso.ux(L,M)+c^2*s^2*(aniso.uxr(L,M)+aniso.uy(L,M))+s^4*aniso.uyr(L,M))*...
        conj(c^4*aniso.uyr(L,M)+c^2*s^2*(aniso.uxr(L,M)+aniso.uy(L,M))+s^4*aniso.ux(L,M)))/2;
    weight(3,L,M) = las*real((c^4*aniso.uy(L,M)+c^2*s^2*(aniso.uyr(L,M)+aniso.ux(L,M))+s^4*aniso.uxr(L,M))*...
        (c^4*aniso.uxr(L,M)+c^2*s^2*(aniso.ux(L,M)+aniso.uyr(L,M))+s^4*aniso.uy(L,M)))/2;
    for M = 2:L
        weight(1,L,M) = real((c^4*aniso.ux(L,M)+c^2*s^2*(aniso.uxr(L,M)+aniso.uy(L,M))+s^4*aniso.uyr(L,M))*...
            (c^4*aniso.uy(L,M)+c^2*s^2*(aniso.ux(L,M)+aniso.uyr(L,M))+s^4*aniso.uxr(L,M)) + ...
            las^2*(c^4*aniso.uxr(L,M)+c^2*s^2*(aniso.uyr(L,M)+aniso.ux(L,M))+s^4*aniso.uy(L,M))*...
            (c^4*aniso.uyr(L,M)+c^2*s^2*(aniso.uxr(L,M)+aniso.uy(L,M))+s^4*aniso.ux(L,M)));
        weight(2,L,M) = las*real((c^4*aniso.ux(L,M)+c^2*s^2*(aniso.uxr(L,M)+aniso.uy(L,M))+s^4*aniso.uyr(L,M))*...
            conj(c^4*aniso.uyr(L,M)+c^2*s^2*(aniso.uxr(L,M)+aniso.uy(L,M))+s^4*aniso.ux(L,M)));
        weight(3,L,M) = las*real((c^4*aniso.uy(L,M)+c^2*s^2*(aniso.uyr(L,M)+aniso.ux(L,M))+s^4*aniso.uxr(L,M))*...
            (c^4*aniso.uxr(L,M)+c^2*s^2*(aniso.ux(L,M)+aniso.uyr(L,M))+s^4*aniso.uy(L,M)));
    end
end
weight = weight/max(weight(:));

if length(p)>5
    epp = zeros(size(t)); epo = epp; eop = epp;
    for L = 1:size(aniso.ux,1)
        for M = 1:L
            epp = epp + weight(1,L,M)*exp(-L*(L-1)*t*rot(1) - (M-1)^2*rot(2));%/(2*L-1);
            epo = epo + weight(2,L,M)*exp(-L*(L-1)*t*rot(1) - (M-1)^2*rot(2));%/(2*L-1);
            eop = eop + weight(3,L,M)*exp(-L*(L-1)*t*rot(1) - (M-1)^2*rot(2));%/(2*L-1);
        end
    end
    zz = [ones(length(t),1) zeros(length(t),4)];
    for j=0:n
%         zz(:,3) = zz(:,3) + epp.*f(t,tau,j*period) + ...
%             (epo.*(t0>=0) + eop.*(t0<0)).*f(t,tau,j*period+delay) + ...
%             (eop.*(t0>=0) + epo.*(t0<0)).*f(t,tau,(j+1)*period-delay);
         zz(:,3) = zz(:,3) + epp.*f(t,tau,j*period);
         zz(:,4) = zz(:,4) + (epo.*(t0>=0) + eop.*(t0<0)).*f(t,tau,j*period+delay);
         zz(:,5) = zz(:,5) + (eop.*(t0>=0) + epo.*(t0<0)).*f(t,tau,(j+1)*period-delay);
    end
else
    zz = [ones(length(t),1) zeros(length(t),1)];
end
for j=-n:n
    zz(:,2) = zz(:,2) + weight(1,1,1)*ff(t0,tau,j*period) + ...
        weight(2,1,1)*ff(t0,tau,j*period+delay) + ...
        weight(3,1,1)*ff(t0,tau,(j+1)*period-delay);
end

if length(p)>5
    c0 = lsqnonneg(zz(:,1:2),x(:,2));
    z0 = zz(:,1:2)*c0;
    %c = lsqnonneg(zz,x(:,1));
    c = zz\x(:,1)
    z = zz*c;
    if nargin>5 && ~isempty(bld)
        plot(t0, x(:,1), 'o', t0, z);
        axis tight
        drawnow
    end
    err = sum(sum((x(:,1)-z).^2))+sum(sum((x(:,2)-z0).^2));
else
    c0 = lsqnonneg(zz(:,1:2),x);
    z0 = zz(:,1:2)*c0;
    c = c0;
    z = z0;
    if nargin>5 && ~isempty(bld)
        plot(t0, x, 'o', t0, z0);
        axis tight
        drawnow
    end
    err = sum(sum((x-z0).^2));
end

% close; para = Simplex('AntibunchDelay',[63 25 16 2 30 2],[0 24 0 0 0 0],[inf 26 inf inf inf 10],[],[],0.032*bin,y(:,10),weight,1);

