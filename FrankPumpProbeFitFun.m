function err = FrankPumpProbeFitFun(param, irf, y2, p, c1, A1, tau1)
%	FrankPumpProbeFitFun(param, irf, y, p) returns the Least-Squares deviation between the data y 
%	and the computed values. 
%
%	  y =  yoffset + A(1)*convol(irf,exp(-t/tau(1)/(1-exp(-p/tau(1)))) + ...
%           A(m+1)*convol(irf,exp(-(t-t1)/tau(m+1)/(1-exp(-p/tau(m+1)))) + ...
%
%	param(1:2) are the color shift values between irf and y.
%	param(3:...) are the decay times.
%	irf is the measured Instrumental Response Function.
%	y is the measured fluorescence decay curve.
%	p is the time between to laser excitations (in number of TCSPC channels).

n = length(irf);
t = 1:n;
tp = (1:p)';
c2 = param(1);
tau = param(2:length(param)); 
tau = tau(:)'; 
tau1 = tau1(:)';
x1 = exp(-(tp-1)*(1./tau1))*diag(1./(1-exp(-p./tau1)));
x2 = exp(-(tp-1)*(1./tau))*diag(1./(1-exp(-p./tau)));
irs1 = (1-c1+floor(c1))*irf(rem(rem(t-floor(c1)-1, n)+n,n)+1) + (c1-floor(c1))*irf(rem(rem(t-ceil(c1)-1, n)+n,n)+1);
irs2 = (1-c2+floor(c2))*irf(rem(rem(t-floor(c2)-1, n)+n,n)+1) + (c2-floor(c2))*irf(rem(rem(t-ceil(c2)-1, n)+n,n)+1);
z1 = convol(irs1, x1)*A1;
z2 = convol(irs2, x2);
z2 = [ones(size(z1,1),1) z1 z2];
A2 = z2\y2;
%A2 = lsqnonneg(z2,y2);
z2 = z2*A2;
if 1
    semilogy(t, irs2/max(irs2)*max(y2), t, y2, 'o', t, z2);
	drawnow
end
err = sum(sum((z2-y2).^2./abs(z2)))/(n-length(tau))


