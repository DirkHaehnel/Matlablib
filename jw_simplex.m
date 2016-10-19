function [x, dx, steps] = jw_simplex(fname, x, xmin, xmax, tol, steps, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10)
%	[x, dx, steps] = SIMPLEX('F', X0, XMIN, XMAX, TOL, STEPS, P1, ..., P10) 
%	attempts to return a vector x and its error dx, so that x minimzes the 
%	function F(x) near the starting vector X0 under the conditions that 
% 	xmin <= x <= xmax.
%	TOL is the relative termination tolerance dF/F; (default = 1e-10)
%	STEPS is the maximum number of steps; (default = 200*number of parameters).
%	The returned value of STEPS is the actual number of performed steps. 
%	SIMPLEX allows for up to 10 additional arguments for the function F.
%	SIMPLEX uses a Nelder-Mead simplex search method.

x = x(:);
if nargin<5
	tol = 1e-10;
if nargin<4
		xmax = Inf*ones(length(x),1);
		if nargin<3
			xmin = -Inf*ones(length(x),1);
		end
	end
elseif isempty(tol)
tol = 1e-5;
end
if nargin<6
	steps = [];
end
if isempty(xmin) xmin = -Inf*ones(size(x)); end
if isempty(xmax) xmax = Inf*ones(size(x)); end
xmin = xmin(:);
xmax = xmax(:);
xmax(xmax<xmin) = xmin(xmax<xmin);
x(x<xmin) = xmin(x<xmin);
x(x>xmax) = xmax(x>xmax);
xfix = zeros(size(x));
xfix(xmin==xmax) = xmin(xmin==xmax);
mask = diag(xmin~=xmax);
mask(:, sum(mask)==0) = [];
tmp = xmin==xmax;
x(tmp) = [];
xmin(tmp) = [];
xmax(tmp) = [];

evalstr = [fname];
evalstr=[evalstr, '(mask*x+xfix'];
for i=1:nargin - 6
	evalstr = [evalstr, ',p', int2str(i)];
end
evalstr = [evalstr, ')'];

n = length(x);
if n==0 
	x = xfix;
	dx = zeros(size(xfix));
	steps = 0;
	%break
    return %jw replace break with return as instructed by matlab error
end
if isempty(steps)
	steps = 200*n;
end

xin = x(:);
v = 0.9*xin;
v(v<xmin) = xmin(v<xmin);
v(v>xmax) = xmax(v>xmax);
x(:) = v; fv = eval(evalstr); 
for j = 1:n
	y = xin;
	if y(j) ~= 0
		y(j) = (1 +.1*rand)*y(j);
	else
		y(j) = 0.1;
	end
	if y(j)>=xmax(j) y(j) = xmax(j); end
	if y(j)<=xmin(j) y(j) = xmin(j); end
	v = [v y];
	x(:) = y; f = eval(evalstr);
	fv = [fv f];
end
[fv, j] = sort(fv);
v = v(:,j);
count = n+1;

% Parameter settings for Nelder-Meade
alpha = 1; beta = 1/2; gamma = 2;

% Begin of Nelder-Meade simplex algorithm
while count < steps
	if 2*abs(fv(n+1)-fv(1))/(abs(fv(1))+abs(fv(n+1))) <= tol
		break
	end

	% Reflection:
	vmean = mean(v(:, 1:n)')';
	vr = (1 + alpha)*vmean - alpha*v(:, n+1);
	x(:) = vr;
	fr = eval(evalstr); 
	count = count + 1; 
	vk = vr; fk = fr;

	%if fr < fv(1) && prod(xmin<=vr) && prod(vr<=xmax) jorg
    if fr < fv(1) && prod(xmin)<=prod(vr) && prod(vr)<=prod(xmax) %jw change, & below
		% Expansion:
		ve = gamma*vr + (1-gamma)*vmean;
		x(:) = ve;
		fe = eval(evalstr);
		count = count + 1;
		if fe < fv(1) && prod(xmin)<=prod(ve) && prod(ve)<=prod(xmax) %change as above
			vk = ve; fk = fe;
		end
	else
		vtmp = v(:,n+1); ftmp = fv(n+1);
		if fr < ftmp && prod(xmin)<=prod(vr) && prod(vr)<=prod(xmax) %change as above
			vtmp = vr; ftmp = fr;
		end
		% Contraction:
		vc = beta*vtmp + (1-beta)*vmean;
		x(:) = vc;
		fc = eval(evalstr); 
		count = count + 1;
		if fc < fv(n) && prod(xmin)<=prod(vc) && prod(vc)<=prod(xmax) %change as above
			vk = vc; fk = fc;
		else
			% Shrinkage:
			for j = 2:n
				v(:, j) = (v(:, 1) + v(:, j))/2;
				x(:) = v(:, j);
				fv(j) = eval(evalstr); 
			end
			count = count + n-1;
			vk = (v(:, 1) + v(:, n+1))/2;
			x(:) = vk;
			fk = eval(evalstr); 
			count = count + 1;
		end
	end
	v(:, n+1) = vk;
	fv(n+1) = fk;
	[fv, j] = sort(fv);
	v = v(:,j);
end

x = v(:,1);
dx = abs(v(:,n+1)-v(:,1));
x = mask*x + xfix;
dx = mask*dx;
if count>=steps
	disp(['Warning: Maximum number of iterations (', int2str(steps),') has been exceeded']);
else
	steps = count;
end
