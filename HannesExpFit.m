function [err, z, c1, c2] = HannesExpFit(p, t, x, kt, intrel, kequi);

keff = p(1);
kp = p(2);
km = kequi*p(2);
kt1 = p(3);

t = t(:); p = p(:)';

M = [-keff-kp, kt, km ,kt1; keff, -kt-kp, 0, 0; kp, 0, -km, 0; 0, kp, 0, -kt1];
ev = eig(M).';
z1 = exp(t*ev);
keff = intrel*keff;
M = [-keff-kp, kt, km ,kt1; keff, -kt-kp, 0, 0; kp, 0, -km, 0; 0, kp, 0, -kt1];
ev = eig(M).';
z2 = exp(t*ev);

c1 = z1\x(:,1);
c2 = z2\x(:,2);
z = [z1*c1 z2*c2];

plot(t, z, t, x, '.'); drawnow

err = abs(sum(sum((x-z).^2)))


