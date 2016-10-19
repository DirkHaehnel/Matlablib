function res = GaussianLaser(r,z,w0)

% Function GaussianLaser(r,z,w0) calculates the intensity distribution of a 
% focused Gaussian laser in water
% the unit of length is chosen in such a way that lambda = 2*pi

nsample = 1.33; % water

R = w0*sqrt(1 + (2*z/w0^2/nsample).^2);

res = w0^2*exp(-2*r.^2./R.^2)./R.^2;