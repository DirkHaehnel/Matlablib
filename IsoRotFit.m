function [err, z] = IsoRotFit(p,theta,y);

% IsoRotFit fits rotational diffusion data for an isotropic rotator

y = y(:);
z = IsoRotFun(theta,p,100)';
c = z\y;
z = c*z;

err = sqrt(sum((y-z).^2))