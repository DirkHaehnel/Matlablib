% Afterpulse Figures

tau = 5;
phi = 0.2;
iquer = 1;
Tv = [1 5.0000001 9];
eps = 0.2;

t = 10.^(0:0.01:2);

T = Tv(1);
acf = iquer^2*(1 + (phi/(1-phi))*exp(-t/tau));
semilogx(t,acf/acf(end),'o','symbolsize',2);
hold on
for j=1:length(Tv)
    T=Tv(j);
    acf = (1+eps)^2*iquer^2 + iquer^2*(phi/(1-phi))*(exp(-t/tau) + eps*(eps+2)/(1-(T/tau)^2)*(exp(-alpha*t)-T/tau*exp(-t/T))) + ...
        eps*iquer/T*exp(-t/T);
    semilogx(t,acf/acf(end),'color',[(j-1)/length(Tv) 0 1-(j-1)/length(Tv)]);
end
hold off
axis tight
xlabel('time [\mus]')
ylabel('autocorrelation')