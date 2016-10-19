function [err z] = CommonFactorLSQ(p,x,y)

for j=1:length(p)
    for k=1:size(y,2)
        x(:,k) = x(:,k) + x(:,k+j*size(y,2))*p(j);
    end
end
x = x(:,1:size(y,2));
for j=1:size(y,2)
    z(:,j) = (x(:,j)\y(:,j))*x(:,j);
end

err = sum(sum(abs((y-z).^2)));
