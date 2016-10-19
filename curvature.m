function [kappa, t, xx, yy] = curvature(x,y)

x=x(:); y=y(:);
n = length(x);
t = (0:n-1)';
m = round(n/2);

xs = polyfit(t,x,m);
ys = polyfit(t,y,m);

xs1 = xs(1:m).*(m:-1:1);
ys1 = ys(1:m).*(m:-1:1);

xs2 = xs1(1:m-1).*(m-1:-1:1);
ys2 = ys1(1:m-1).*(m-1:-1:1);

t = 0:0.1:n-1;
kappa = [];
tmp1 = conv(xs1,ys2)-conv(xs2,ys1);
tmp2 = conv(xs1,xs1) + conv(ys1,ys1);
kappa = [kappa polyval(tmp1,t)./polyval(tmp2,t).^(3/2)];
xx = polyval(xs,t);
yy = polyval(ys,t);
t=t/(n-1);



