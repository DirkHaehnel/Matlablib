function [err,z,zz,M] = ExonucleaseFun(p, t, seq, y, ind)

n = length(seq);
if nargin<5 || isempty(ind)
    ind = 1:n; ind = ind(seq==max(seq));
end

if length(p(:))<max(seq)^2
    tmp = p; p = zeros(max(seq)); 
    for j=1:max(seq)
        p(j,1:max(seq)) = tmp(j); 
    end
else
    p = reshape(p,max(seq),max(seq));
end
p    
for j=1:length(ind) weight(j) = p(seq(ind(j)),seq(ind(j)+1)); end
z0 = zeros(n, 1);
z0(1) = 1;
M = zeros(n,n);
for j=1:n
    if j<n
        M(j,j) = -p(seq(j),seq(j+1));
        M(j+1,j) = p(seq(j),seq(j+1));
    end
end
for k=1:length(t)
    z(:,k) = expm(t(k)*M)*z0;
end;

if nargin>3 && ~isempty(y)
    zz = (weight*z(ind,:))';
    c = lsqnonneg(zz,y(:));
    zz = zz*c;
    plot(t,y,'o',t,zz); drawnow
    err = sum((y(:)-zz).^2./(zz+(zz==0)));
else
    err = z;
    z = M;
    zz = weight;
end

return

% model data

t = 0:2e-3:1;
k0 = 108; sig = 66;
kv = 10:300;
weight = exp(-(kv-k0).^2/2/sig^2);
seq = [1, 2, 2, 1, 4, 1, 2, 2, 1, 1, 2, 1, 2, 1, 2, 1, 2, 1, 1, 1, 2, 2, 2, 2, 1, 1, 2, 2, 1, 2, 1, 2, 2, 2, 2, 1, 1, 4, 3];

z = zeros(length(seq),length(t));
for j=1:length(kv)
    z = z + weight(j)*ExonucleaseFun(kv(j), t, ones(size(seq)));
end

seq = [1, 2, 2, 1, 4, 1, 2, 2, 1, 1, 2, 1, 2, 1, 2, 1, 2, 1, 1, 1, 2, 2, 2, 2, 1, 1, 2, 2, 1, 2, 1, 2, 2, 2, 2, 1, 1, 4, 3];