t = juliandate(2010,1:61,1);
plot(t,1.04.^((t-t(1))/365.2422)*350)
set(gca,'XTick',t(1:12:end))
set(gca,'XTickLabel',2010:2015)
xlabel('Jahresbeginn')
ylabel('Schuldlast (k€)')
axis tight
grid
tt = interp1(1.04.^((t-t(1))/365.2422)*350,t,400,'spline');
hold on; plot(tt,1.04.^((tt-t(1))/365.2425)*350,'ob'); hold off

[y,m,d] = jd2date(tt);