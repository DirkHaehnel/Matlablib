function err = fluofitfit(param, irf, y, p)
%	FluoFitFun(param, irf, y, p) returns the Least-Squares deviation between the data y 
%	and the computed values. 
%	FluoFitFun assumes a function of the form:
%
%	  y =  offset + A(1)*convol(irf,exp(-t/tau(1)/(1-exp(-p/tau(1)))) + ...
%
%	param(1) is the color shift value between irf and y.
%	param(2) is the offset value.
%	param(3:...) are the decay times.
%	irf is the measured Instrumental Response Function.
%	y is the measured fluorescence decay curve.
%	p is the time between to laser excitations (in number of TCSPC channels).


global plotdata
n = length(irf);
t = 1:n;
tp = 1:p;
c = param(1);
offset = param(2);
tau = param(3:length(param));
x = diag(1./(1-exp(-p./tau)))*exp(-(1./tau)*(tp-1));
x = x';
z = zeros(n, length(tau));
irs = (1-c+floor(c))*irf(rem(rem(t-floor(c)-1, n)+n,n)+1) + (c-floor(c))*irf(rem(rem(t-ceil(c)-1, n)+n,n)+1);
for j=1:length(tau)
	z(:,j) = convol(irs, x(:,j))';
end
A = abs(z\(y-offset));
z = z*A + offset;
set(plotdata(end),'ydata',z);
drawnow
err = (z-y)'*((z-y)./abs(z))/(n-length(tau));
