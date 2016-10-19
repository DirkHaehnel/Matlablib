function [err, c, z] = PaidFun(p, x, y)

p = reshape(p,2,length(p)/2); row = ones(1,size(p,2)); 
x = x(:); col = ones(size(x));
y = y(:);

z = sqrt(p(2,:)/2).*(1+p(1,:)./p(2,:));
z = mOrtho(4,max(x),z);
z = z(:,x+1)'; 
z = exp(x/2*log(p(2,:)/2)-col*(p(1,:)+p(2,:)/2)).*z./(gamma(x+1)*row); 
z(z<1e-10 | isnan(z))=0;
%c = lsqnonneg(z,y);
%c = z(z>0)\y(z>0);
c = z\y;
loglog(x,y,'o',x,z*c); drawnow
err = sum((y-z*c).^2)

