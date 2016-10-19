function y = IsoRotFun(theta,t,maxn);

% IsoRotFun computes the angle-change distribution for an isotropic rotator

if nargin<3
    maxn = 5;
end

y = zeros(1,length(theta));
for j=0:maxn
    tmp = legendre(j,cos(theta));
    y = y + sqrt((2*j+1)/4/pi)*tmp(1,:)*exp(-j*(j+1)*t);
end



