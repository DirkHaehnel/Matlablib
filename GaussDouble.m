function [err, c, z, A] = GaussDouble(p,t,y,bck);

if isempty(t)
	t = 1:length(y);
end
t = t(:);
y = y(:);
nn = length(p)/3;

A = zeros(length(t),2*nn);
for j = 1:nn
   A(:,2*j-1) = exp(-(t-p(end/3+j)-p(j)).^2/p(2*end/3+j)^2/2);
   A(:,2*j) = exp(-(t-p(j)+p(end/3+j)).^2/p(2*end/3+j)^2/2);
end
if nargin==4
    A(:,2*nn+1) = bck(:);
end
ind = y>0;
% c = lsqnonneg(A(ind,:),y(ind));
c = A(ind,:)\y(ind);
z = A*c;
plot(t,y,'o',t,z); drawnow
err = sum((z(ind)-y(ind)).^2./abs(z(ind)));

