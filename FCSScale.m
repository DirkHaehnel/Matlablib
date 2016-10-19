function [y, z] = FCSScale(x, k, j)

if nargin<2
    k=10;
end
if nargin<3
    j=10;
end

y = (x(k:end,:)-ones(size(x(k:end,1)))*x(end,:))./(ones(size(x(k:end,1)))*(abs(mean(x(k:k+j,:))-x(end,:))));

if nargout>1
	z = (x-ones(size(x,1),1)*x(end,:))./(ones(size(x,1),1)*(abs(mean(x(k:k+j,:))-x(end,:))));
end
