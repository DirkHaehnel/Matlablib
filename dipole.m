function [v,pc,ps] = Dipole(theta,z,n)

% [v,pc,ps] = Dipole(theta,z,n) calculates the dipole radiation along emission angle theta of a dipole at distance z 
% from an interface with refractive index n = [n(1) n(2)]

z = z(:)';
col = ones(size(z));
theta = abs(theta(:));
ind = theta<=pi/2;

v = zeros(length(theta),length(z));
pc = v; ps = v;

tmp = theta(ind);
w1 = n(1)*cos(tmp);
w2 = sqrt(n(2)^2-n(1)^2*sin(tmp).^2);
rp = (w1*n(2)^2-w2*n(1)^2)./(w1*n(2)^2+w2*n(1)^2);
rs = (w1-w2)./(w1+w2);
vtmp = n(1)*(sin(tmp).^2*col).*abs(1+(rp*col).*exp(2*i*w1*z)).^2;
pctmp = (w1.^2/n(1)*col).*abs(1+(rp*col).*exp(2*i*w1*z)).^2;
pstmp = n(1)*abs(1+(rs*col).*exp(2*i*w1*z)).^2;
v(ind) = vtmp;
pc(ind) = pctmp;
ps(ind) = pstmp;
tmp = theta(~ind);
w1 = -sqrt(n(1)^2-n(2)^2*sin(tmp).^2);
w2 = n(2)*cos(tmp);
tp = 2*n(1)*n(2)*w1./(w1*n(2)^2+w2*n(1)^2);
ts = 2*w1./(w1+w2);
vtmp = n(2)^3/n(1)^2*((w2.^2./abs(w1).^2.*sin(tmp).^2.*abs(tp).^2)*col).*exp(2*imag(w1)*z);
pctmp = n(2)^3/n(1)^4*((w2.^2.*abs(tp).^2)*col).*exp(2*imag(w1)*z);
pstmp = n(2)^3/n(1)^2*((w2.^2./abs(w1).^2.*abs(ts).^2)*col).*exp(2*imag(w1)*z);
v(~ind) = vtmp;
pc(~ind) = pctmp;
ps(~ind) = pstmp;
