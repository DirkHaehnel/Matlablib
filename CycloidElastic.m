y0 = 0.1:0.1:1;
tmax = 20;

opts = odeset('RelTol',1e-12,'AbsTol',1e-20);

for j=1:length(y0)
    [t,y] = ode15s('CycloidElasticFun',[0 20],[y0(j) 0]',opts);
    plot(t,y); drawnow
    Tper(j) = mean(diff(t(maxl(y(:,1)))));
end



