function [rp,rs,tp,ts] = DielMirror(theta,ni,nf,n,a,nn)

% [rp,rs,tp,ts] = DielMirror(theta,ni,nf,n,a,nn) calculates the reflection and transmission 
% of an dielectric mirror with alternatng refractive [n(1) n(2)], thickness [a(1) a(2)] and nn layers 
% ni - refractive index of incident medium
% nf - refractive index of substrate

wi = ni*cos(theta);
w1 = sqrt(n(1)^2-ni^2*sin(theta).^2);
w2 = sqrt(n(2)^2-ni^2*sin(theta).^2);
wf = sqrt(nf^2-ni^2*sin(theta).^2);

Mi11 = ni/2./wi;
Mi12 = 1/2/ni*ones(size(wi));
Mi21 = -Mi11;
Mi22 = Mi12;

Mf1 = wf/nf;
Mf2 = nf;

M11 = cos(a(1)*w1).*cos(a(2)*w2)-(n(2)^2*w1)./(n(1)^2*w2).*sin(a(1)*w1).*sin(a(2)*w2);
M22 = cos(a(1)*w1).*cos(a(2)*w2)-(n(1)^2*w2)./(n(2)^2*w1).*sin(a(1)*w1).*sin(a(2)*w2);
M12 = -i*(w1/n(1)^2.*cos(a(2)*w2).*sin(a(1)*w1)+w2/n(2)^2.*cos(a(1)*w1).*sin(a(2)*w2));
M21 = -i*(n(1)^2./w1.*cos(a(2)*w2).*sin(a(1)*w1)+n(2)^2./w2.*cos(a(1)*w1).*sin(a(2)*w2));

dM = sqrt(4*M12.*M21+(M11-M22).^2);
l1 = (M11 + M22 + dM)/2;
l2 = (M11 + M22 - dM)/2;

eig11 = M12;
eig21 = l1-M11;
eig12 = l2-M22;
eig22 = M21;

tmp = eig11.*eig22-eig12.*eig21;
inv11 = eig22./tmp.*l1.^nn;
inv12 = -eig12./tmp.*l1.^nn;
inv21 = -eig21./tmp.*l2.^nn;
inv22 = eig11./tmp.*l2.^nn;

r1 = inv11.*Mf1 + inv12.*Mf2;
r2 = inv21.*Mf1 + inv22.*Mf2;
tmp = eig11.*r1+eig12.*r2;
r2 = eig21.*r1+eig22.*r2;
r1 = tmp;
rp = (Mi21.*r1+Mi22.*r2)./(Mi11.*r1+Mi12.*r2);
tp = 1./(Mi11.*r1+Mi12.*r2);


Mi11 = 1/2*ones(size(wi));
Mi12 = 1/2./wi;
Mi21 = Mi11;
Mi22 = -Mi12;

Mf1 = ones(size(wi));
Mf2 = wf;

M11 = cos(a(1)*w1).*cos(a(2)*w2)-w2./w1.*sin(a(1)*w1).*sin(a(2)*w2);
M22 = cos(a(1)*w1).*cos(a(2)*w2)-w1./w2.*sin(a(1)*w1).*sin(a(2)*w2);
M12 = -i*(1./w1.*cos(a(2)*w2).*sin(a(1)*w1)+1./w2.*cos(a(1)*w1).*sin(a(2)*w2));
M21 = -i*(w1.*cos(a(2)*w2).*sin(a(1)*w1)+w2.*cos(a(1)*w1).*sin(a(2)*w2));

dM = sqrt(4*M12.*M21+(M11-M22).^2);
l1 = (M11 + M22 + dM)/2;
l2 = (M11 + M22 - dM)/2;

eig11 = M12;
eig21 = l1-M11;
eig12 = l2-M22;
eig22 = M21;

tmp = eig11.*eig22-eig12.*eig21;
inv11 = eig22./tmp.*l1.^nn;
inv12 = -eig12./tmp.*l1.^nn;
inv21 = -eig21./tmp.*l2.^nn;
inv22 = eig11./tmp.*l2.^nn;

r1 = inv11.*Mf1 + inv12.*Mf2;
r2 = inv21.*Mf1 + inv22.*Mf2;
tmp = eig11.*r1+eig12.*r2;
r2 = eig21.*r1+eig22.*r2;
r1 = tmp;
rs = (Mi21.*r1+Mi22.*r2)./(Mi11.*r1+Mi12.*r2);
ts = 1./(Mi11.*r1+Mi12.*r2);

