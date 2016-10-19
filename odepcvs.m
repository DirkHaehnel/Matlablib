function [tout, yout] = odepcvs(F, t0, tfinal, y0, tol, trace, k)
% ODEPCVS  Variable-Step Predictor-Corrector Ordinary Differential Equation 
%        Solver
%
% [tout, yout] = odepcvs(F, t0, tfinal, y0, tol, trace, k)
%
% INPUTS
% F     - String containing name of user-supplied problem description.
%         Call: yprime = fun(t,y) where F = 'fun'.
%         t      - Time (scalar).
%         y      - Solution column-vector.
%         yprime - Returned derivative column-vector; yprime(i) = dy(i)/dt.
% y0     - Initial value column vector
% t0     - Initial value of t
% tfinal - Final value of t
% tol    - Relative tolerance for Local Truncation Error (default = 1.e-6)
% trace - If nonzero, each step is printed. (Default: trace = 0).
% k      - Order of predictor-corrector method to be used (default = 4)
%
% OUTPUTS
% tout   - Returned integration time points (row-vector)
% yout   - Returned solution, one solution column-vector per tout-value
%
% The resulting trajectory can be displayed by: plot(tout, yout).
% Written by Heath Hofmann, October 3rd, 1996
% heath@eecs.berkeley.edu

if nargin < 7, k = 4; end
if nargin < 6, trace = 0; end
if nargin < 5, tol = 1e-6; end

m = length(y0);

%Determine coefficents
gamma = zeros(k+1,1);
gammas = zeros(k+1,1);
gamma(1) = 1;
for i = 2:k+1, for j = 1:i-1, gamma(i) = gamma(i) - gamma(j)/(i+1-j); end, end
for i=1:k+1, for j = 1:i, gammas(i) = gammas(i) + gamma(j); end,end
beta = (pascal(k,1)'*gamma(1:k))';
betak = beta(1);
beta(1:k-1) = beta(2:k);
beta(k) = 0;
betas = (pascal(k,1)'*gammas(1:k))';

deltas = (betas - beta)/betak;
cpp1 = gamma(k+1);
cspp1 = gammas(k+1);
W = cpp1/(cspp1-cpp1);
pow = 1/8;

% Initialization
hmax = (tfinal - t0);
hmin = (tfinal - t0)/1e6;
h = max(tol/100, hmin);
y = y0;
yout = y0;
t = t0;tout = t;
f = [];

PI = kron(abs(pascal(k+1,1))',eye(m));
Qb = zeros(k+1);
Qb(1,1)=1;
Qb(2,2)=1;
for i = 2:k+1, for j=3:k+1, Qb(j,i) = (i-1)*(2-j)^(i-2); end, end
Qb = inv(Qb);

Q = kron(Qb, eye(m));
G = kron([betak 1 zeros(1,k-1)]',eye(m));
Gt = Q*G;

if trace
   clc, t, y
end

%First-order startup procedure

withintol=1;
count = 0;

while (t < tfinal) & (h >= hmin) & (count < k) 

   %First-order Adams Bashforth (Euler's method), Predictor
   fp = feval(F,t,y)
   yp = y + h*fp;   

   %First-order Adams-Moulton, Corrector
   fc = feval(F,t,yp);
   yc = y + h*fc;

   %Determine error using "Milne's trick"
   gamma1 = -.5*(yp-yc);

   % Estimate the error and the acceptable error
   delta = norm(gamma1,'inf');
   tau = tol*max(norm(y,'inf'),1.0);
   withintol = (delta <= tau);
   % Update the solution only if the error is acceptable
   if withintol
      t = t + h;
      y = yc;
      tout = [tout t];
      yout = [yout y];
      count = count+1;
      f = [fc f];
   else
      count = 0;
      tout = t0;
      yout = y0;
      h = min(hmax, 0.8*h*(tau/delta)^pow);
   end
   if trace
      home, t, y
   end
end

d = f*deltas';
Y = [y; h*d];
for i = 1:k-1, Y = [Y; h*f(:,i)]; end
Zc = Q*Y;

% Predictor-Corrector Method

while (h >= hmin) & (t < tfinal)

   %Predict
   Zp = PI*Zc;

   %Correct
   fnpk = feval(F,t,Zp(1:m));
   Fn = (h*fnpk - Zp(m+1:2*m));
   Zc = Zp + Gt*Fn;

   %Determine error using Milne's trick

   error = W*(Zc(1:m) - Zp(1:m));

   % Estimate the error and the acceptable error
   delta = norm(error,'inf');
   tau = tol*max(norm(y,'inf'),1.0);

   % Update the solution only if the error is acceptable
   if delta <= tau
      t = t + h;
      y = Zc(1:m);
      tout = [tout t];
      yout = [yout y];
   end
   if trace
      home, t, y
   end

   % Update the step size
   if delta ~= 0.0
      hold = h;
      h = min(hmax, 0.8*h*(tau/delta)^pow);
      if t + h > tfinal, h = tfinal - t; end
      alpha = h/hold;
   else
      alpha = 1;
   end

   %Adjust step size
   for i = 1:k
      Zc((i*m+1):(i+1)*m) = alpha^i*Zc((i*m+1):(i+1)*m);
   end

end

if (t < tfinal) 
   disp('SINGULARITY LIKELY.')
   t
end

