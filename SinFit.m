function [err,c,z] = SinFit(p,t,x)

t = t(:); x = x(:);
z = sin(2*pi/p(1)*(t-p(2)));
c = [ones(length(x),1) t z]\x;
z = [ones(length(x),1) t z]*c;

plot(t,x,'o',t,z); drawnow

err = sum(abs(x-z).^2)


