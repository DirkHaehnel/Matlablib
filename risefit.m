function [err, y, t] = risefit(p,x)

x = x(:)';
t = (0:length(x)-1);
y = exp(p(1)*t./(p(2)+t.^p(3)));
%y = atan(p(1)*(t-p(2).^p(3)));
y = [y; ones(1,length(x))];

c = x/y;
y = c*y;

plot(t,x,'o',t,y); drawnow
err = sum((x-y).^2)