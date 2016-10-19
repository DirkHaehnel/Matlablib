function [tau, cx] = DistFit(x,dt);

close all
if nargin<2
    dt = 0.25;
end

N = 100;

x(1) = [];
t = dt*(1:length(x))';
tau = 5./exp((0:N)/N*log(max(t))); % distribution of decay times
M = [ones(size(t)) exp(-t*tau)]; 

cx = lsqnonneg(M,x);
z = M*cx;
for j=1:size(M,2)
   M(:,j) = M(:,j)./sqrt(z);
end
cx = lsqnonneg(M,sqrt(x));
z = (M*cx).^2;
semilogy(t,x,'.b',t,z,'g','linewidth',2); times16

v = axis;
v(1) = min(t);
v(2) = max(t);
axis(v);
xlabel('time in \mus');
ylabel('log count');
figure;
subplot(2,1,1);
plot(t,(x-z)./sqrt(z)); times16
v = axis;
v(1) = min(t);
v(2) = max(t);
axis(v);
xlabel('time in \mus');
ylabel('weighted residual');

ind=1:length(cx)-2;
len = length(ind);
tau = 1./tau;
fac = sqrt(tau(1:end-1)/tau(2:end));
subplot(2,1,2)
semilogx(reshape([fac*tau(ind);fac*tau(ind);tau(ind)/fac;tau(ind)/fac],4*len,1),reshape([0*tau(ind);cx(ind+1)';cx(ind+1)';0*tau(ind)],4*len,1));
patch(reshape([fac*tau(ind);fac*tau(ind);tau(ind)/fac;tau(ind)/fac],4*len,1),reshape([0*tau(ind);cx(ind+1)';cx(ind+1)';0*tau(ind)],4*len,1),'b');
times16
xlabel('decay time in \mus');
ylabel('distribution');