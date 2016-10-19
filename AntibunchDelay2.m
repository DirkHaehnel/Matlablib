function [err, c1, c2, z, zz] = AntibunchDelay2(p, t0, x, aniso, las, bld)

% p(1) ... null time
% p(2) ... period
% p(3) ... delay
% p(4) ... decay time
% p(5:6) ... rotation time

t0 = t0(:)-p(1);
period = p(2);
t = abs(t0);
delay = mod(p(3),period);
tau = p(4);
rot = p(5);
del = p(6);
if nargin<5 || isempty(las) 
    las = 1;
end
n = 2*ceil(max(t)/period)*ceil(tau/period);

f = inline('sinh(x/tau).*exp(-delay/tau).*(x<=delay) + sinh(delay/tau).*exp(-x/tau).*(x>delay)','x','tau','delay');
ff = inline('exp((x-delay)/tau).*(x<=delay)/2 + exp(-(x-delay)/tau).*(x>delay)/2','x','tau','delay');

weight(1,:) = real([aniso.ux00*conj(aniso.uy00)/2 aniso.ux20*conj(aniso.uy20)/2, -aniso.ux21*conj(aniso.uy21) aniso.ux22*conj(aniso.uy22) ...
    aniso.ux40*conj(aniso.uy40)/2, -aniso.ux41*conj(aniso.uy41) aniso.ux42*conj(aniso.uy42), -aniso.ux43*conj(aniso.uy43) aniso.ux44*conj(aniso.uy44)]);

weight(2,:) = real([aniso.ux00*conj(aniso.uy00r)/2 aniso.ux20*conj(aniso.uy20r)/2, -aniso.ux21*conj(aniso.uy21r) aniso.ux22*conj(aniso.uy22r) ...
    aniso.ux40*conj(aniso.uy40r)/2, -aniso.ux41*conj(aniso.uy41r) aniso.ux42*conj(aniso.uy42r), -aniso.ux43*conj(aniso.uy43r) aniso.ux44*conj(aniso.uy44r)]);

weight(3,:) = real([aniso.ux00r*conj(aniso.uy00)/2 aniso.ux20r*conj(aniso.uy20)/2, -aniso.ux21r*conj(aniso.uy21) aniso.ux22r*conj(aniso.uy22) ...
    aniso.ux40r*conj(aniso.uy40)/2, -aniso.ux41r*conj(aniso.uy41) aniso.ux42r*conj(aniso.uy42), -aniso.ux43r*conj(aniso.uy43) aniso.ux44r*conj(aniso.uy44)]);

weight(4,:) = (real([aniso.ux00*conj(aniso.ux00)/2 aniso.ux20*conj(aniso.ux20)/2, -aniso.ux21*conj(aniso.ux21) aniso.ux22*conj(aniso.ux22) ...
    aniso.ux40*conj(aniso.ux40)/2, -aniso.ux41*conj(aniso.ux41) aniso.ux42*conj(aniso.ux42) -aniso.ux43*conj(aniso.ux43) aniso.ux44*conj(aniso.ux44)]) + ...
    real([aniso.uy00*conj(aniso.uy00)/2 aniso.uy20*conj(aniso.uy20)/2, -aniso.uy21*conj(aniso.uy21) aniso.uy22*conj(aniso.uy22) ...
    aniso.uy40*conj(aniso.uy40)/2, -aniso.uy41*conj(aniso.uy41) aniso.uy42*conj(aniso.uy42), -aniso.uy43*conj(aniso.uy43) aniso.uy44*conj(aniso.uy44)]))/2;

weight(5,:) = real([aniso.ux00*conj(aniso.ux00r)/2 aniso.ux20*conj(aniso.ux20r)/2, -aniso.ux21*conj(aniso.ux21r) aniso.ux22*conj(aniso.ux22r) ...
    aniso.ux40*conj(aniso.ux40r)/2, -aniso.ux41*conj(aniso.ux41r) aniso.ux42*conj(aniso.ux42r), -aniso.ux43*conj(aniso.ux43r) aniso.ux44*conj(aniso.ux44r)]);

weight = diag([1+las^2 las las 1+las^2 las])*weight;

epp = weight(1,1) + weight(1,2)*exp(-t/rot) +  weight(1,3)*exp(-t/rot-t/del) + weight(1,4)*exp(-t/rot-4*t/del) + ...
    weight(1,5)*exp(-10/3*t/rot) + weight(1,6)*exp(-10/3*t/rot-t/del) + weight(1,7)*exp(-10/3*t/rot-4*t/del) + ...
    weight(1,8)*exp(-10/3*t/rot-9*t/del) + weight(1,9)*exp(-10/3*t/rot-16*t/del);

epo = weight(2,1) + weight(2,2)*exp(-t/rot) +  weight(2,3)*exp(-t/rot-t/del) + weight(2,4)*exp(-t/rot-4*t/del) + ...
    weight(2,5)*exp(-10/3*t/rot) + weight(2,6)*exp(-10/3*t/rot-t/del) + weight(2,7)*exp(-10/3*t/rot-4*t/del) + ...
    weight(2,8)*exp(-10/3*t/rot-9*t/del) + weight(2,9)*exp(-10/3*t/rot-16*t/del);

eop = weight(3,1) + weight(3,2)*exp(-t/rot) +  weight(3,3)*exp(-t/rot-t/del) + weight(3,4)*exp(-t/rot-4*t/del) + ...
    weight(3,5)*exp(-10/3*t/rot) + weight(3,6)*exp(-10/3*t/rot-t/del) + weight(3,7)*exp(-10/3*t/rot-4*t/del) + ...
    weight(3,8)*exp(-10/3*t/rot-9*t/del) + weight(3,9)*exp(-10/3*t/rot-16*t/del);

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
    %plot(t0, z, 'b', t0, x, 'or'); 
    plot(t0, x, 'o', t0, z); 
    drawnow
end

err = sum(sum((x-z).^2)./z./t);
%err = sum(sum((x-z).^2));

% close; para = Simplex('AntibunchDelay',[63 25 16 2 30 2],[0 24 0 0 0 0],[inf 26 inf inf inf 10],[],[],0.032*bin,y(:,10),weight,1);

