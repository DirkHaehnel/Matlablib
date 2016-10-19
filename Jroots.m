function y=jroots(m, n);
% Function y = jroots(m, n) finds the first n roots of the function BesselJ(m, x)

y(1) = abs(simplex('jfun', 0, [], [], [], [], m));
if y(1)==0
	[tmp, xl] = max(besselj(m, 1:2*m)); xl=xl/pi;
	[tmp, xr] = min(besselj(m, 1:2*m)); xr=xr/pi;
	y(2) = simplex('jfun', (xr+xl)/2, xl, xr, [], [], m)
	for k=3:n y(k)=simplex('jfun', y(k-1)+1, y(k-1)+0.5, y(k-1)+1.5, [], [], m); end
else
	for k=2:n y(k)=simplex('jfun', y(k-1)+1, y(k-1)+0.5, y(k-1)+1.5, [], [], m); end
end
