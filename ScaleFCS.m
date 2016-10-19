function y = ScaleFCS(x, k, j);

if nargin<2
    k=10;
end
if nargin<3
    j=10;
end

y = (x(k:end,:)-mean(x(end,:)))./(ones(size(x,1)-k+1,1)*(mean(x(k:k+j,:))-mean(x(end,:))));

