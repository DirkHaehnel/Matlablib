function [c, offset, A, tau, dc, doffset, dtau, irs, z, chi] = fluofit(irf, y, p, dt, m, options, fitfun, pic)
% The function FLUOFIT performs a fit of a multi-exponential decay curve.
% It is called by: 
% [c, off, A, tau, dc, doffset, dtau, irs, z, err] = fluofit(irf, y, p, dt, m, opt, fun, pic).
% The function arguments are:
% irf 	= 	Instrumental Response Function
% y 	= 	Fluorescence decay data
% p 	= 	Time between laser exciation pulses (in nanoseconds)
% dt 	= 	Time width of one TCSPC channel (in nanoseconds)
% m 	= 	Number of desired decay components
% opt	=	If opt(1)==1, then colorshift is included, else not.
%		If opt(2)==1, then offset is included, else not.
%		If length(opt)>2, then opt(3:...) are used as initial guess times, 
%		and no automatic initial guessing is made.
% fun	=	If fun=='lsfit', then a least squares fit is applied; 
%		if fun=='mlfit', then a maximum likelihood fit is applied.
% pic	=	If pic==1, than figures are drawn, else not.
%
% The return parameters are:
% c	=	Color Shift (time shift of the IRF with respect to the fluorescence curve)
% off	=	Offset
% A	=	Amplitudes of the different decay components
% tau	=	Decay times of the different decay components
% err	=	Final value of the fit function
% dc	=	Color shift error
% doffset	= 	Offset error
% dtau	=	Decay times error
% irs	=	IRF, shifted by the value of the colorshift
% z	=	Fitted fluorecence curve
% 
% The program needs the following m-files: simplex.m, lsfit.m, mlfit.m, and convol.m.
% (c) 1996 Jörg Enderlein

global bild
if nargin==8
	if pic>0 bild = 1; else bild = 0; end
else
	bild = 1;
end
if nargin>=7
	if isempty(fitfun) fitfun = 'lsfit'; end
else
	fitfun = 'lsfit';
end
close all
p = p/dt;
irf = irf(:)';
y = y(:)';
n = length(irf); 
tp = 1:p;
t = 1:length(y);
init = 1;
shifton = 1;
offseton = 1;
if nargin>5
	tmp = length(options);
	if tmp>0
	if options(1)==0  shifton = options(1); end
	end
	if tmp>1
		if options(2)==0  offseton = options(2); end
	end
	if tmp>2
		init = 0;
		tau = options(3:tmp);
	end
end
if shifton
	sh_min = -5;
	sh_max = 5;
else
	sh_min = 0;
	sh_max = 0;
end
if init
	param(1) = m;
	for i=2:m 
		param(i) = -param(i-1)*(m-i+1)/i;
	end;
	B = y(m+1:n)';
	for i=1:m
		B = B + param(i)^2*y(m+1-i:n-i)';
	end
	for i=1:m
		tmp = -param(i)*y(m+1:n-i)';
		for j=1:m-i
			tmp = tmp + param(j)*param(j+i)*y(m+1-j:n-j-i)';
		end
		B = [[tmp; zeros(i,1)] B [zeros(i,1); tmp]];
	end
	d = -m:m;
	V = spdiags(B, d, n-m, n-m);
	D = [y(m:n-1)'];
	for i=2:m
		D = [D y(m+1-i:n-i)'];
	end
	D = D + 1e-10*(D==0);
	D(:,m+1:2*m+1)=ones(n-m,m+1);
	yy = y(m+1:n)';
	for k=sh_min:sh_max
		irs=irf(rem(rem(t-k-1, n)+n,n)+1);
		for i=1:m
			D(:,m+i) = irs(m+2-i:n+1-i)';
		end
		VD = V\D;
		param = (VD'*D)\(VD'*yy);
		tau = roots([param(m:-1:1)' -1]);
		if ~isempty(tau)
			for i=1:m
				tmp = (1-tau(i)./tau);
				tmp(i) = [];
				A(i) = polyval(param(2*m:-1:m+1), tau(i))/prod(tmp);
			end
			tau = 1./log(tau);
			offset = param(2*m+1)/(1-sum(param(1:m)));
			z = A*diag(1./(1-exp(-p./tau)))*exp(-(1./tau)*(tp-1));
			z = convol(irs, z)' + offset;
            size(y)
            size(z)
			tmp = ((y-z)./abs(z))*(y-z)'/(n-m);
			chi(k-sh_min+1) = tmp;
		else
			chi(k-sh_min+1) = nan;
		end
	end
	[chi, c] = min(chi);
	c = c + sh_min - 1;
	irs=irf(rem(rem(t-c-1, n)+n,n)+1);
	for i=1:m
		D(:,m+i) = irs(m+2-i:n+1-i)';
	end
	VD = V\D;
	param = (VD'*D)\(VD'*yy);
	tau = roots([param(m:-1:1)' -1]);
	for i=1:m
		tmp = (1-tau(i)./tau);
		tmp(i) = [];
		A(i) = polyval(param(2*m:-1:m+1), tau(i))/prod(tmp);
	end
	tau = abs(1./log(tau))
	offset = param(2*m+1)/(1-sum(param(1:m)));
	z = A*diag(1./(1-exp(-p./tau)))*exp(-(1./tau)*(tp-1));
	z = convol(irs, z) + offset;
else
	tau = tau(:)
	m = length(tau);
	x = diag(1./(1-exp(-p./tau)))*exp(-(1./tau)*(tp-1));
	z = zeros(m, n);
	for j=1:m
		z(j,:) = convol(irf, x(j,:));
	end
	A = y/z;
	z = A*z;
	c = 0;
	offset = options(2);
end
if offseton offset=abs(offset); else offset = 0; end
disp([c, offset, A, tau']);
if bild
	hold on
	plot(t, irf/max(irf)*max(y), t, y, 'bo' , 'EraseMode', 'none');
	global Plothandle
	Plothandle = plot(t,z,'EraseMode','xor');
	drawnow
end
disp('Fit =                Parameters =');
param = [c; offset; tau];
% Decay times and Offset are assumed to be positive.
paramin = [-Inf, 0, zeros(1,length(tau))];
paramax = [Inf, Inf, Inf*ones(1,length(tau))];
if ~shifton  paramin(1) = 0; paramax(1) = 0; end
if ~offseton  paramax(2) = 0; end
[param, dparam] = simplex(fitfun, param, paramin, paramax, 1e-8, [], irf(:), y(:), p);
c = param(1);
dc = dparam(1);
offset = param(2);
doffset = dparam(2);
tau = param(3:length(param));
dtau = dparam(3:length(param));
x = diag(1./(1-exp(-p./tau)))*exp(-(1./tau)*(tp-1));
irs = (1-c+floor(c))*irf(rem(rem(t-floor(c)-1, n)+n,n)+1) + (c-floor(c))*irf(rem(rem(t-ceil(c)-1, n)+n,n)+1);
for j =1:m
	z(j,:) = convol(irs, x(j,:));
end
A = (z*z')\(z*(y'-offset));
z = A'*z + offset;
chi = ((y-z)./abs(z))*(y-z)'/(n-m);
t = dt*t;
tau = dt*tau;
c = dt*c;
dtau = dt*dtau;
dc = dt*dc;
disp([chi A' tau']);
if bild
	hold off
	plot(t,log10(y),t,log10(irf),t,log10(z));
	v = axis;
	v(1) = min(t);
	v(2) = max(t);
	axis(v);
	xlabel('Time in ns');
	ylabel('Log Count');
	s = sprintf('COF = %3.3f   %3.3f', c, offset);
	text(max(t)/2,v(4)-0.05*(v(4)-v(3)),s);
	s = ['AMP = '];
	A = A/sum(A);
	for i=1:length(A)
		s = [s sprintf('%1.3f',A(i)) '   '];
	end
	text(max(t)/2,v(4)-0.12*(v(4)-v(3)),s);
	s = ['TAU = '];
	for i=1:length(tau)
		s = [s sprintf('%3.3f',tau(i)) '   '];
	end
	text(max(t)/2,v(4)-0.19*(v(4)-v(3)),s);
	figure;
	subplot(2,1,1);
	plot(t,(y-z)./sqrt(abs(z)));
	v = axis;
	v(1) = min(t);
	v(2) = max(t);
	axis(v);
	if fitfun=='lsfit'
		s = sprintf('CHI2 = %3.2f', chi);
	else
		s = sprintf('MLE = %3.2f', chi);
	end
	xlabel('Time in ns');
	ylabel('Residue');
	title(s);
end
