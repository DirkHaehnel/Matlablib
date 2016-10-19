function [corr, t] = RotoDiffPotential(t,drot,eps,al,nmax)
global x al n m k j drot eps

% drot = [dpara, dortho]

if nargin<4 || isempty(al)
    al = 0;
end

if nargin<5 || isempty(nmax)
    nmax = 180;
end

x = linspace(0,pi,nmax);
M = RotoDiffPotentialMat(acos(x),al);
corr = ones(length(t),4);
for n = 1:1%4 % order of correlator
    for j = 0:n % element of vector of n-th order correlator
        for m = -n:n % final z-component of momentum
            for k = 0:n; % initial z-component of momentum
                if ~all(M{n}(:,m+n+1,k+1,j+1)<1e-15)
                    [n m k]
                    u = RotoDiffPotentialPDE(x,t);
                    corr(:,n) = corr(:,n) + binomial(n,j)*(u*M{n}(:,m+n+1,k+1,j+1));
                    if k>0
                        corr(:,n) = corr(:,n) + binomial(n,j)*(u*M{n}(:,-m+n+1,k+1,j+1));
                    end
                end
            end
        end
    end
end
for j=1:4
    corr(:,j) = corr(:,j)/corr(1,j);
end


% --------------------------------------------------------------
function u = RotoDiffPotentialPDE(x,t)
% global x t al pp qq m k j drot
sol = pdepe(0,@rdppde,@rdpic,@rdpbc,x,t);
u = sol(:,:,1);
% --------------------------------------------------------------
function [c,f,s] = rdppde(x,t,u,DuDx)
global m k drot eps
c = sin(x);
f = sin(x).*(DuDx + u*eps*(1-cos(x))); 
s = -u.*(drot(2)*k^2 + 0.5*sum(drot)*m^2 - 2*drot(2)*m*k*cos(x) - 0.5*diff(drot)*m^2*cos(2*x))./sin(x);
% --------------------------------------------------------------
function u0 = rdpic(x)
global eps n j al
if eps>0
    u0 = eps*exp(-eps*(1-cos(x)))/(1-exp(-2*eps));
else
    u0 = ones(size(x))/2;
end
u0 = sin(x).*u0.*sin(x-al).^j.*cos(x-al).^(n-j);
% --------------------------------------------------------------
function [pl,ql,pr,qr] = rdpbc(xl,ul,xr,ur,t)
pl = 0;
ql = 1;
pr = 0;
qr = 1;

