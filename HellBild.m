NA = 1.4;
n = 1.51;
over = [NA*3e3 inf 5e3];
lamex = 0.514;
[fx0, fx2, fz, rho, z, exc, excx, excz] = NewExc([-0.005 0.8], 0, NA, n, n, n, [], 0, [], lamex, over, 0, [], lamex/0.01);
exc = exc(:,1)';
rho = rho(:,1)';
clear fx0 fx2 fz z excx ecxz

exc0 = 1e7;
tau = 1e-9;
kisc = 1e8;
kph = 1e5;
t = (0:0.1:10)*1e-6;
ind = [1 5 9 17 33];
rhov = rho(ind);
zerfall = Hell(exc0*exc(ind),t,tau,kisc,kph);
plot(1e6*t,zerfall'./(ones(size(zerfall,2),1)*zerfall(:,1)'))
xlabel('time [\mus]')
ylabel('rel. amplitude')

tv = [0 0.1 0.2 0.4 0.8 1.6 3.2 6.4]*1e-6;
profile = Hell(exc0*[fliplr(exc) exc],tv,tau,kisc,kph);
plot([-fliplr(rho) rho],profile/exc0)

tv = [0 0.1 0.2 0.4 0.8 1.6 3.2 6.4]*1e-6;
profile = Hell(exc0*exc,tv,tau,kisc,kph);
plot(rho,profile/exc0)
hold on
for j=1:length(rhov)
    plot([rhov(j) rhov(j)],[0 1],':');
end 
hold off
xlabel('\itx\rm [\mum]')
ylabel('rel. amplitude')
axis([0 0.5 0 1])
set(gca,'yaxislocation','right')