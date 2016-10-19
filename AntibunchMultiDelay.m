function [err, c, z] = AntibunchMultiDelay(p, t0, x, bld)

% p(1) ... null time
% p(2) ... period
% p(3) ... delay
% p(4) ... decay time
% p(5) ... rotation time
% p(6) ... average numer of molecules in focus
% p(7) ... delay between measurements
% p(8:9) ... rel brightness values

t0 = t0(:)-p(1);
t = abs(t0);
n = ceil(max(t)/p);
period = p(2);
delay = p(3);
tau = p(4);
rot = p(5);
smd = p(6);
if length(p)>6
    delta = p(7);
else
    delta = 0;
end


f = inline('sinh(x/tau).*exp(-delay/tau).*(x<delay) + sinh(delay/tau).*exp(-x/tau).*(x>=delay)','x','tau','delay');
ff = inline('exp(-abs(x-delay)/tau)','x','tau','delay');

col = ones(length(x),1);
para = 4+5*exp(-t/rot);
ortho = 4-exp(-t/rot);

kk = 9;
delay = mod(delay,period);
for k=1:size(x,2)
    zpp = smd*para.*f(t,tau,0) + smd*(smd-1)*ff(t,tau,0);
    zpo = smd*ortho.*f(t,tau,delay) + smd*(smd-1)*ff(t,tau,delay);
    zop = smd*ortho.*f(t,tau,period-delay) + smd*(smd-1)*ff(t,tau,period-delay);
    for j=1:n
        zpp = zpp + smd*para.*f(t,tau,j*period) + smd*(smd-1)*ff(t,tau,j*period);
        zpo = zpo + smd*ortho.*f(t,tau,j*period+delay) + smd*(smd-1)*ff(t,tau,j*period+delay);
        zop = zop + smd*ortho.*f(t,tau,(j+1)*period-delay) + smd*(smd-1)*ff(t,tau,(j+1)*period-delay);
    end
    tmp = zpo;
    zpo(t0<0) = zop(t0<0);
    zop(t0<0) = tmp(t0<0);
    delay = mod(delay + delta,period);
    %zz = zpp + p(8)*zpo + p(9)*zop;
    %zz = 6*zpp + 9*zpo + zop;
    %zz = 6*zpp + zpo + 9*zop;
    %zz = zpp + zpo + zop;
    %c(:,k) = lsqnonneg([ones(size(x,1),1) zz],x(:,k));
    %c(:,k) = [ones(size(x,1),1) zz]\x(:,k);
    %z(:,k) = [ones(size(x,1),1) zz]*c(:,k);
    c(:,k) = lsqnonneg([ones(size(x,1),1) zpp zpo zop],x(:,k));
    z(:,k) = [ones(size(x,1),1) zpp zpo zop]*c(:,k);
    if k==kk
        %tst = [z(:,k) c(2,k)*zpp c(2,k)*p(8)*zpo c(2,k)*p(9)*zop;];
        tst = [z(:,k) c(2,k)*zpp c(2,k)*zpo c(2,k)*zop;];
    end
end
    
if nargin>3 && ~isempty(bld)
    plot(t0, x(:,kk), 'o', t0, z(:,kk), 'markersize', 2.5); drawnow
    %plot(t0, x(:,kk), 'o', t0, tst, 'markersize', 2.5); drawnow
end

err = sum(sum((x-z).^2./abs(z)));



