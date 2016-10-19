function [v,pc,ps] = DipoleE(theta,z,n)

% [v,pc,ps] = DipoleE(theta,z,n) calculates the electric field amplitudes of dipole radiation along emission angle theta 
% of a dipole at distance z from an interface with refractive index n = [n(1) n(2)]

z = z(:)';
col = ones(size(z));
theta = abs(theta(:));
ind = theta<pi/2;

v = zeros(length(theta),length(z));
pc = v; ps = v;

if sum(ind)>0
    tmp = theta(ind);
    w1 = sqrt(n(1)^2-n(2)^2*sin(tmp).^2);
    w2 = n(2)*cos(tmp);
    tp = 2*n(1)*n(2)*w1./(w1*n(2)^2+w2*n(1)^2);
    ts = 2*w1./(w1+w2);
    fac = sqrt(n(2))*n(2)*w2./w1;
    v(ind,:) = ((fac.*sin(tmp).*n(2)/n(1).*tp)*col).*exp(i*w1*z);
    pc(ind,:) = ((fac.*w1/n(1).*tp)*col).*exp(i*w1*z);
    ps(ind,:) = ((fac.*ts)*col).*exp(i*w1*z);
end
