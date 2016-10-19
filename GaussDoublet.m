function [err, c1, c2, z1, z2, A1, A2] = GaussDoublet(p,t,y,pic)

polyorder = 1;

[n, m] = size(y);
if m>n 
    y = y';
    [n, m] = size(y);
end

if isempty(t)
	t = 1:n;
end
t = t(:);
nn = length(p)/3;

A1 = zeros(length(t),nn + polyorder);
A2 = A1;
for j = 1:nn
   A1(:,j) = exp(-(t-p(j)).^2/p(2*end/3+j)^2/2);
   A2(:,j) = exp(-(t-p(end/3+j)).^2/p(2*end/3+j)^2/2);
end
for j=1:polyorder
    A1(:,nn+j) = t.^(j-1);
    A2(:,nn+j) = t.^(j-1);
end
ind = sum(y,2)>0 & all(isfinite(y),2);
c1 = lsqnonneg(A1(ind,:),y(ind,1));
c2 = lsqnonneg(A1(ind,:),y(ind,2));
z1 = A1*c1;
z2 = A2*c2;
if nargin>3 && ~isempty(pic)
    plot(t,y,'o',t,z1,t,z2); drawnow
end
%err = sum((z1(ind)-y(ind,1)).^2./abs(z1(ind))) + sum((z2(ind)-y(ind,2)).^2./abs(z2(ind)));
err = sum((z1(ind)-y(ind,1)).^2.) + sum((z2(ind)-y(ind,2)).^2.);

