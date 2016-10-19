function [c, offset, A, tau, irs, z, zz, t, chi] = jw_tcspcfit(irf, y, p, dt, tau)
%2007-04-27: John W modified Jorg's code so it would work in MATLAB 2007a
%changes are noted in comment lines below
%
%The function TCSPCFIT performs a fit of a multi-exponential decay curve.
% It is called by: 
% [c, offset, A, tau, dc, doffset, dtau, irs, z, t, chi] = fluofit(irf, y, p, dt, tau, init).
% The function arguments are:
% irf 	= 	Instrumental Response Function
% y 	= 	Fluorescence decay data
% p 	= 	Time between laser exciation pulses (in nanoseconds)
% dt 	= 	Time width of one TCSPC channel (in nanoseconds)
% tau 	= 	Initial guess times
%
% The return parameters are:
% c	=	Color Shift (time shift of the IRF with respect to the fluorescence curve)
% offset	=	Offset
% A	    =   Amplitudes of the different decay components
% tau	=	Decay times of the different decay components
% dc	=	Color shift error
% doffset	= 	Offset error
% dtau	=	Decay times error
% irs	=	IRF, shifted by the value of the colorshift
% z	    =	Fitted fluorecence curve
% t     =   time axis
% chi   =   chi2 value
% 
% The program needs the following m-files: simplex.m, lsfit.m, mlfit.m, and convol.m.
% (c) 2004 Jörg Enderlein

%fitfun is passed to jorg's simplex function where it is assembled
%into a string calling the lsfit function
%jorg's code 'irflsfit' calls a non-existent function 'irflsfit'
%fitfun = 'irflsfit'; jorg's code
fitfun = 'lsfit'; %jw change because there is no fxn 'irflsfit', only 'lsfit'

close all                   
irf = irf(:);               
offset = 0;                  
y = y(:);                   
n = length(irf);                       
c = 0;
p = p/dt;                   
tp = (1:p)';                
tau = tau(:)'/dt;            
t = 1:length(y);            
m = length(tau);            
x = exp(-(tp-1)*(1./tau))*diag(1./(1-exp(-p./tau)));
irs = (1-c+floor(c))*irf(rem(rem(t-floor(c)-1, n)+n,n)+1) + (c-floor(c))*irf(rem(rem(t-ceil(c)-1, n)+n,n)+1);
z = convol(irs, x);
z=z(:);
z = [ones(size(z,1),1) irs z];
A = z\y;
z = z*A;
close all

hold on
plot(t, irf/max(irf)*max(y), t, y, 'bo' , 'EraseMode', 'none');
global Plothandle
Plothandle = plot(t,z,'EraseMode','xor');
drawnow
disp('Fit =                Parameters =');
param = [c; tau'];
% Decay times and Offset are assumed to be positive.
paramin = [-Inf, zeros(1,length(tau))];
paramax = [Inf, Inf*ones(1,length(tau))];
% **************************************************************
[param, dparam, steps] = jw_simplex(fitfun, param, paramin, paramax, [], [], irf(:), y(:), p);
% **************************************************************
c = param(1);
dc = dparam(1);
tau = param(2:length(param))';
dtau = dparam(2:length(param));
x = exp(-(tp-1)*(1./tau))*diag(1./(1-exp(-p./tau)));
irs = (1-c+floor(c))*irf(rem(rem(t-floor(c)-1, n)+n,n)+1) + (c-floor(c))*irf(rem(rem(t-ceil(c)-1, n)+n,n)+1);
z = convol(irs, x);
z=z(:);
z = [ones(size(z,1),1) irs z];
z = z./(ones(n,1)*sum(z));
A = z\y;
for j=1:length(A)
    zz(:,j) = z(:,j)*A(j);
end
z = z*A;
dtau = dtau;
dc = dt*dc;
chi = sum((y-z).^2./abs(z))/(n-m);
t = dt*t;
tau = dt*tau';
c = dt*c;
offset = A(1); A(1) = [];
hold off
subplot('position',[0.1 0.4 0.8 0.5])
plot(t,log10(y),t,log10(irf),t,log10(z)); %blue  green   red
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
subplot('position',[0.1 0.1 0.8 0.2])
plot(t,(y-z)./sqrt(abs(z)));
v = axis;
v(1) = min(t);
v(2) = max(t);
axis(v);
xlabel('Time in ns');
ylabel('Residue');
s = sprintf('%3.3f', chi);
text(max(t)/2,v(4)-0.1*(v(4)-v(3)),['\chi^2 = ' s]);
set(gcf,'units','normalized','position',[0.01 0.05 0.98 0.83])
