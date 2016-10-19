function res = D3CEF(r,z,a,NA,nsample)

% Function D3CEF(r,z,a,NA,nsample) calculates the collection efficiency
% the unit of length is chosen in such a way that lambda = 2*pi

if size(r,1)==1 | size(r,2)==1
    [z,r]=meshgrid(z,r);
end

w0 = 2/tan(asin(NA/nsample))/nsample

R = w0*sqrt(1 + (2*pi*z/pi/w0^2/nsample).^2);

res = zeros(size(r));

t = logical(r<=abs(a-R));
res(t) = pi*R(t).^2;
res(t & a<R) = pi*a^2;
t = logical(abs(a-R)<r & r<R+a);
theta1 = acos((a^2+r(t).^2-R(t).^2)./r(t)/a/2);
theta2 = acos((R(t).^2+r(t).^2-a^2)./r(t)./R(t)/2);
res(t) = theta1*a^2 + theta2.*R(t).^2 - sqrt((a+R(t)+r(t)).*(-a+R(t)+r(t)).*(a-R(t)+r(t)).*(a+R(t)-r(t)))/2;

res = res./R.^2/pi*a^2/max(a,w0)^2;