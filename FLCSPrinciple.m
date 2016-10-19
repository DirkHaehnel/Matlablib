% FLCS principle

tau = [1 2];
t = (-2:0.01:12)';

irf = exp(-7*t.^2);
f1 = 0*irf; f2 = f1;
for j=1:length(irf)
    f1(j:end) = f1(j:end) + irf(j)*exp(-t(j:end)/tau(1));
    f2(j:end) = f2(j:end) + irf(j)*exp(-t(j:end)/tau(2));
end
f1 = f1/sum(f1);
f2 = f2/sum(f2);

subplot(3,3,2)
plot(t,f1/max(f1),'r')
axis tight
subplot(3,3,3)
plot(t,f2/max(f2),'b')
axis tight

ff = inv([f1 f2]'*[f1 f2])*[f1 f2]';

subplot(3,3,4)
plot(t,ff(1,:),'r')
axis tight
subplot(3,3,7)
plot(t,ff(2,:),'b')
axis tight
