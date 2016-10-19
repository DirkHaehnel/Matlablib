function [x,flag,relres,iter,resvec] = pcgnull(A,tol,maxit,M1,M2,x0)
%PCGNULL    Preconditioned Conjugate Gradients Method
%   X = PCGNULL(A) attempts to solve the homogenious system of linear equations A*X=0
%   for X such that NORM(X)=1.  
%   The coefficient matrix A must be Hermitian and positive
%   semi-definite. If the null-space of A is more than one dimentional,
%   i.e. A*X=0 allows multiple linear independent solutions,
%   PCGNULL will still converge to one solution, namely, the orthoprojection
%   of the initial guess to the null-space.
%   If B is a symmetric matrix with known minimal eigenvalue lambda,
%   X = PCGNULL(B - lambda .* speye(size(B)) attemps to compute
%   an eigenvector of B, corresponding to lambda.
%   PCGNULL will start iterating from an initial
%   guess which by default is a random vector of length N.  Iterates
%   are produced until the method either converges, fails, or has
%   computed the maximum number of iterations.  Convergence is achieved
%   when an iterate X has relative residual NORM(A*X)/NORM(X) less than
%   or equal to the tolerance of the method.  The default tolerance is
%   1e-6.  The default maximum number of iterations is the minimum of
%   N and 20.  No preconditioning is used.
%
%   PCGNULL(A,TOL) specifies the tolerance of the method.
%
%   PCGNULL(A,TOL,MAXIT) specifies the maximum number of iterations.
%
%   PCGNULL(A,TOL,MAXIT,M) and PCG(A,TOL,MAXIT,M1,M2) use Hermitian
%   positive definite (left) preconditioner M or M=M1*M2 and effectively
%   solve the system inv(M)*A*X = 0 for X.  If any preconditioner is
%   given as the empty matrix ([]) then that preconditioner is considered to
%   be the identity matrix.  Since systems of the form M*Y = R are solved for Y
%   using \ (backslash), it is wise to factor a non-triangular preconditioner
%   M into its Cholesky factors i.e. replace  PCGNULL(A,TOL,MAXIT,M)  with
%   M2 = chol(M);  M1 = M2';  PCGNULL(A,TOL,MAXIT,M1,M2).
%
%   PCGNULL(A,TOL,MAXIT,M1,M2,X0) specifies the initial guess.  If X0 is given
%   as the empty matrix ([]) then the default random vector is used.
%
%   X = PCGNULL(A,TOL,MAXIT,M1,M2,X0) returns a solution X.  If PCG
%   converged, then a message to that effect is displayed.  If PCG
%   failed to converge after the maximum number of iterations or halted
%   for any reason, then a warning message is printed displaying the
%   relative residual NORM(A*X)/NORM(X) and the iteration number at
%   which the method stopped or failed.
%
%   [X,FLAG] = PCGNULL(A,TOL,MAXIT,M1,M2,X0) returns a solution X and
%   a FLAG which describes the convergence of PCG.  If FLAG is
%    0 then PCG converged to the desired tolerance TOL within MAXIT
%      iterations without failing for any reason.
%    1 then PCG iterated MAXIT times but did not converge.
%    2 then a system of equations of the form M*Y = R was ill-conditioned
%      and did not return a useable result when solved by \ (backslash).
%    3 then PCG stagnated (two consecutive iterates were the same).
%    4 then one of the scalar quantities calculated during PCG became
%      too small or too large to continue computing.
%   Whenever FLAG is not 0, the solution X returned is that with
%   minimal norm residual computed over all the iterations.
%   No messages are displayed if the FLAG output is specified.
%
%   [X,FLAG,RELRES] = PCGNULL(A,TOL,MAXIT,M1,M2,X0) also returns the
%   relative residual NORM(A*X)/NORM(X).  If FLAG is 0, RELRES <= TOL.
%
%   [X,FLAG,RELRES,ITER] = PCGNULL(A,TOL,MAXIT,M1,M2,X0) also returns the
%   iteration number at which X was computed: 0 <= ITER <= MAXIT.
%
%   [X,FLAG,RELRES,ITER,RESVEC] = PCGNULL(A,TOL,MAXIT,M1,M2,X0) also returns a
%   vector of the relative residual norms at each iteration, starting from
%   RESVEC(1) = NORM(A*X0)/NORM(X0).  If FLAG is 0, then RESVEC is of length ITER+1
%   and RESVEC(END) <= TOL.
%
%   Example:
%     A = delsq(numgrid('L',25));
%     R = cholinc(A,1e-3);
%     shift = min(eig(A));
%     A = A - shift .* speye(size(A)); %this produces a symmetric positive
%     %semi-definite matrix A
%     x0 = ones(length(A),1);
%     [x,flag,relres,iter,resvec] = pcgnull(A,1e-8,70,[],[],x0);
%   flag is 1 since PCG will not converge to within the requested tolerance of
%   1e-8 (confirmed by the value of relres) in 70 iterations
%   without preconditioning.
%     [x,flag,relres,iter,resvec] = pcgnull(A,1e-8,70,R',R,x0);
%   flag is 0 since PCG will converge to within the requested tolerance of
%   1e-8 (confirmed by the value of relres) in 12 iterations
%   when preconditioned by the
%   incomplete Cholesky factorization with a drop tolerance of 1e-3.  You may
%   follow the progress of PCG by plotting the relative residuals at each
%   iteration starting from the initial guess (iterate number 0) with
%     semilogy(0:iter,resvec,'-o')
%
%   See also PCG, BICG, BICGSTAB, CGS, GMRES, QMR, SLASH, CHOLINC.

% Revision 1.0, 1999, by Andrew Knyazev, andrew.knyazev@na-net.ornl.gov
% Tested under MATLAB version 5.1-5.3.  
% This is a modified version of PCG, Revision 1.6, 1996, by Penny Anderson, 
% modified with the permission of The MathWorks, Inc., the copyright owner.  
% This Revision 1.0 may not be used with any products other than products of 
% The MathWorks, Inc., nor may it be used in or as part of another computer program. 
% MATLAB is a registered trademark of The MathWorks, Inc.

% Check for an acceptable number of input arguments
if nargin < 1
  es = sprintf(['Not enough input arguments.']);
  error(es);
end
if nargin > 6
  es = sprintf(['Too many input arguments.']);
  error(es);
end

% Check matrix inputs have appropriate sizes
[m,n] = size(A);
if m ~= n
  error('Matrix must be square.');
end

% Assign default values to unspecified parameters
if nargin < 2
  tol = 1e-6;
end
if nargin < 3
  maxit = min(n,20);
end
if nargin >= 4 & ~isempty(M1)
  [mM1,nM1] = size(M1);
  if ~(mM1==m & nM1==m)
    es = sprintf(['Preconditioner must be a square matrix' ...
                  ' of size %d to match the coefficient matrix.'],m);
    error(es);
  else
    existM1 = 1;
  end
else
  existM1 = 0;
end
if nargin >= 5 & ~isempty(M2)
  [mM2,nM2] = size(M2);
  if ~(mM2==n & nM2==n)
    es = sprintf(['Preconditioner must be a square matrix' ...
                  ' of size %d to match the coefficient matrix.'],n);
    error(es);
  else
    existM2 = 1;
  end
else
  existM2 = 0;
end
if nargin == 6 & ~isempty(x0)
  [mx0,nx0] = size(x0);
  if ~(mx0==n & nx0==1)
    es = sprintf(['Initial guess must be a column vector of' ...
                  ' length %d to match the coefficient matrix.'],n);
    error(es);
  else
    if ~(norm(x0) == 0)
      es = sprintf(['Warning: Initial guess must be a nonzero column vector.' ...
                  ' A random initial guess will now be used.'],n);
      x = randn(n,1);
    else
    x = x0;
    end
  end
else
  x = randn(n,1);  
end

% Set up for the method
flag = 1;
xmin = x;                          % Iterate which has minimal residual so far
imin = 0;                          % Iteration at which xmin was computed
tolb = tol;                  % Relative tolerance
r = - A * x;                     % Zero-th residual
normr = norm(r) / norm (x);         % Relative Norm of residual

if normr <= tolb                   % Initial guess is a good enough solution
  flag = 0;
  relres = normr;
  iter = 0;
  resvec = normr;
  if nargout < 2
    os = sprintf(['The initial guess has relative residual %0.2g' ...
                  ' which is within\nthe desired tolerance %0.2g' ...
                  ' so PCG returned it without iterating.'], ...
                  relres,tol);
    disp(os);
  end
  return;
end

resvec = zeros(maxit+1,1);         % Preallocate vector for norm of residuals
resvec(1) = normr;                 % resvec(1) = norm(A*x0)/norm(x0)
normrmin = normr;                  % Norm of minimum residual
rho = 1;
stag = 0;                          % stagnation of the method

% loop over maxit iterations (unless convergence or failure)

for i = 1 : maxit
  if existM1
%    y = M1 \ r;
y = M1 * r;
%AK
   if isinf(norm(y,inf))
      flag = 2;
      break;
    end 
  else
    y = r;
  end
  if existM2
    z = M2 \ y;
    if isinf(norm(z,inf))
      flag = 2;
      break;
    end
  else
    z = y;
  end
  rho1 = rho;
  rho = r' * z;
  if rho == 0 | rho == Inf
    flag = 4;
    break;
  end
  if i == 1
    p = z;
  else
    beta = rho / rho1;
    if beta == 0 | beta == Inf
      flag = 4;
      break;
    end
    p = z + beta * p;
  end
  q = A * p;
  pq = p' * q;
  if pq == 0 | pq == Inf
    flag = 4;
    break;
  else
    alpha = rho / pq;
  end
  if alpha == Inf
    flag = 4;
    break;
  end
  if alpha == 0                    % stagnation of the method
    stag = 1;
  end

% Check for stagnation of the method
  if stag == 0
    stagtest = zeros(n,1);
    ind = (x ~= 0);
    stagtest(ind) = p(ind) ./ x(ind);
    stagtest(~ind & p ~= 0) = Inf;
    if abs(alpha)*norm(stagtest,inf) < eps
      stag = 1;
    end
  end

  x = x + alpha * p;               % form new iterate
  normr = norm(A*x)/norm(x);
  resvec(i+1) = normr;

  if normr <= tolb                 % check for convergence
    flag = 0;
    iter = i;
    break;
  end

  if stag == 1
    flag = 3;
    break;
  end

  if normr < normrmin              % update minimal norm quantities
    normrmin = normr;
    xmin = x;
    imin = i;
  end

  r = r - alpha * q;

end                                % for i = 1 : maxit

% returned solution is first with minimal residual
if flag == 0
  relres = normr;
else
  x = xmin;
  iter = imin;
  relres = normrmin;
end

x = x ./ norm(x);

% truncate the zeros from resvec
if flag <= 1 | flag == 3
  resvec = resvec(1:i+1);
else
  resvec = resvec(1:i);
end

% only display a message if the output flag is not used
if nargout < 2
  switch(flag)
    case 0,
      os = sprintf(['PCG converged at iteration %d to a' ...
                    ' solution with relative residual %0.2g'], ...
                    iter,relres);

    case 1,
      os = sprintf(['PCG stopped after the maximum %d iterations' ...
                    ' without converging to the desired tolerance %0.2g' ...
                    '\n         The iterate returned (number %d)' ...
                    ' has relative residual %0.2g'], ...
                    maxit,tol,iter,relres);

    case 2,
      os = sprintf(['PCG stopped at iteration %d' ...
                    ' without converging to the desired tolerance %0.2g' ...
                    '\n         because the system involving the' ...
                    ' preconditioner was ill conditioned.' ...
                    '\n         The iterate returned (number %d)' ...
                    ' has relative residual %0.2g'], ...
                    i,tol,iter,relres);

    case 3,
      os = sprintf(['PCG stopped at iteration %d' ...
                    ' without converging to the\n         desired' ...
                    ' tolerance %0.2g because the method stagnated.' ...
                    '\n         The iterate returned (number %d)' ...
                    ' has relative residual %0.2g'], ...
                    i,tol,iter,relres);

    case 4,
      os = sprintf(['PCG stopped at iteration %d' ...
                    ' without converging to' ...
                    ' the desired tolerance %0.2g' ...
                    '\n         because a scalar quantity became too' ...
                    ' small or too large to continue computing.' ...
                    '\n         The iterate returned (number %d)' ...
                    ' has relative residual %0.2g'], ...
                    i,tol,iter,relres);

  end                              % switch(flag)
  if flag == 0
    disp(os);
  else
    warning(os);
  end
end                                % if nargout < 2
