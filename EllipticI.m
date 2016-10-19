function [F,E,Z] = elliptici(u,m,tol)
%   [F,E,Z] = elliptici(u,m,tol) where u is a phase in radians, 0<m<1 is 
%   the modul and tol is the tolerance (optional). Delault value for 
%   the tolerance is eps = 2.220e-16.
%
%   Elliptici uses the method of the arithmetic-geometric mean described 
%   in [1] to determine the value of the Incomplete Elliptic Integrals 
%   of the First, Second Kind and Jacobi's Zeta Function.
%
%       F(u,m) = int(1/sqrt(1-m*sin(t)^2), t=0..u);
%       E(u,m) = int(sqrt(1-m*sin(t)^2), t=0..u);
%       Z(u,m) = E(u,m) - E(m)/K(m)*F(u,m).
%
%   References:
%   [1] M. Abramowitz and I.A. Stegun, "Handbook of Mathematical
%       Functions" Dover Publications", 1965, Ch. 16-17.6.
%   [2] D. F. Lawden, "Elliptic Functions and Applications"
%       Springer-Verlag, vol. 80, 1989
%
%   For support reply: moiseev[at]sissa.it

if nargin<3, tol = eps; end
if nargin<2, error('Not enough input arguments.'); end

if ~isreal(u) | ~isreal(m)
    error('Input arguments must be real.')
end

[mm,nm] = size(m);
[mu,nu] = size(u);
if length(m)==1, m = m(ones(size(u))); end
if length(u)==1, u = u(ones(size(m))); end
if ~isequal(size(m),size(u)), error('U and M must be the same size.'); end

mmax = prod(size(u));

F = zeros(size(u));
E = F;
Z = E;
m = m(:).';    % make a row vector
u = u(:).';

if any(m < 0) | any(m > 1), 
  error('M must be in the range 0 <= M <= 1.');
end

% pre-allocate space and augment if needed
chunk = 10;
a = zeros(chunk,mmax);
c = a;
b = a;
a(1,:) = ones(1,mmax);
c(1,:) = sqrt(m);
b(1,:) = sqrt(1-m);
n = zeros(1,mmax);
i = 1;
while any(abs(c(i,:)) > tol)
    i = i + 1;
    if i > size(a,1)
      a = [a; zeros(chunk,mmax)];
      b = [b; zeros(chunk,mmax)];
      c = [c; zeros(chunk,mmax)];
    end
    a(i,:) = 0.5 * (a(i-1,:) + b(i-1,:));
    b(i,:) = sqrt(a(i-1,:) .* b(i-1,:));
    c(i,:) = 0.5 * (a(i-1,:) - b(i-1,:));
    in = find((abs(c(i,:)) <= tol) & (abs(c(i-1,:)) > tol));
    if ~isempty(in)
      [mi,ni] = size(in);
      n(in) = ones(mi,ni)*(i-1);
    end
end
phin = zeros(i,mmax);
e = zeros(max(n),length(n));
phin(1,:) = u;
i = 0;
while i < n
    i = i + 1;
    in = find(n >= i);
    phin(i+1,:) = phin(i,:);
    if ~isempty(in)
        phin(i+1,in) = atan(b(i,in)/a(i,in)*tan(phin(i,in))) + ...
            pi.*floor((phin(i,in)+0.5*pi)./pi) + phin(i,in);
    end
    e(i,in) = 2.^(i-1);
end
Kk = pi/2./a(i,:);                                                          % Full Eliptic Integral of the First Kind
Ee = Kk .* (1 - 1/2*sum(e(1:n,:).*c(1:n,:).^2,1));                          % Full Eliptic Integral of the Second Kind

F(:) = phin(i,:) ./ (a(i,:)*2^(i-1));                                       % Incomplete Ell. Int. of the First Kind
Z(:) = sum(c(2:n,:).*sin(phin(2:n,:)),1);                                   % Jacobi Zeta Function
E(:) = Z + Ee./Kk.*F;;                                                      % Incomplete Ell. Int. of the Second Kind


