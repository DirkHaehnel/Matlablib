% translational and rotational diffusion of ellipsoids
% MuellerRotoDiff.pdf

close all

qq = 0.5:0.01:2;
dpara = zeros(size(q)); dortho = dpara;
ind = qq<1;
q = qq(ind);
dpara(ind) = 3*q.^2./2./(1-q.^4).*((2-q.^2)./sqrt(1-q.^2).*log((1+sqrt(1-q.^2))./q)-1);
dortho(ind) = 3./4./(1-q.^4).*(q.^2.*(1-2*q.^2)./sqrt(1-q.^2).*log((1+sqrt(1-q.^2))./q)+1);
ind = qq>1;
q = qq(ind);
dpara(ind) = 3*q.^2./2./(q.^4-1).*((q.^2-2)./sqrt(q.^2-1).*atan(sqrt(q.^2-1))+1);
dortho(ind) = 3./4./(q.^4-1).*(q.^2.*(2.*q.^2-1)./sqrt(q.^2-1).*atan(sqrt(q.^2-1))-1);
ind = qq==1;
dpara(ind) = 1;
dortho(ind) = 1;

q = 1:0.01:2;

ft1=sqrt(1-1./q.^2).*q.^(2/3)./log((1+sqrt(1-1./q.^2)).*q); ft2=sqrt(q.^2-1)./q.^(2/3)./atan(sqrt(q.^2-1));

% plot(q,ft1,1./q,ft2)

% prolate elliposid
fra1 = 4*(1-1./q.^2)./(3*(2-2./q.^(4/3)./ft1)); % symmetry axis
frb1 = 4*(1-1./q.^4).*q.^2./(3*(2*q.^(2/3).*(2-1./q.^2)./ft1-2)); % perpendicular axis
% oblate elliposid
fra2 = 4*(1-q.^2)./(3*(2-2*q.^(4/3)./ft2));
frb2 = 4*(1-q.^4)./q.^2./(3*(2./q.^(2/3).*(2-q.^2)./ft2-2));

% figure
% plot(q,4*(1-1./q.^2)./(3*(2-2./q.^(4/3)./ft1)),q,4*(1-1./q.^4).*q.^2./(3*(2*q.^(2/3).*(2-1./q.^2)./ft1-2)),1./q,4*(1-q.^2)./(3*(2-2*q.^(4/3)./ft2)),1./q,4*(1-q.^4)./q.^2./(3*(2./q.^(2/3).*(2-q.^2)./ft2-2))); 
% grid

% rotational diffusion coefficients
figure
plot(q, 1./fra1, 'r', q, 1./frb1, 'b', 1./q, 1./fra2, 'r', 1./q, 1./frb2, 'b');
xlabel('aspect ratio (symmetry axis/transverse axis)')
ylabel('rotational diffusion coefficient (\itD\rm_0)')
legend({'symmetry axis','transverse axis'},2)
grid

% mean diffusion coefficient
figure
plot(q, 1./ft1,'r', q, 1/3*(2./frb1+1./fra1),'b', 1./q, 1./ft2, 'r', 1./q, (2./frb2+1./fra2)/3,'b');
xlabel('aspect ratio (symmetry axis/transverse axis)')
ylabel('mean diffusion coefficient (\itD\rm_0)')
legend({'translational diffusion','rotational diffusion'},4)
grid

% apparent hydro rad trans + rot
figure
plot(q, ft1,'r', q, (3./(2./frb1+1./fra1)).^(1/3), 'b', 1./q, ft2, 'r', 1./q, (3./(2./frb2+1./fra2)).^(1/3),'b');
xlabel('aspect ratio (symmetry axis/transverse axis)')
ylabel('apparent hydrodynamic radius ( \itR\rm_0)')
legend({'translational diffusion','rotational diffusion'},1)
grid

% apparent hydro rad rot
figure
plot(q, (3./(2./frb1+1./fra1)).^(1/3), 'b', 1./q, (3./(2./frb2+1./fra2)).^(1/3),'b');
xlabel('aspect ratio (symmetry axis/transverse axis)')
ylabel('apparent hydrodynamic radius ( \itR\rm_0)')
grid

% ratio of apparent hydro rad
plot(q, ft1./(3./(2./frb1+1./fra1)).^(1/3), 'r', 1./q, ft2./(3./(2./frb2+1./fra2)).^(1/3),'r');
xlabel('aspect ratio (symmetry axis/transverse axis)')
ylabel('\itR_T\rm /\itR_R\rm')
grid


% figure
% plot(q,3*ft1./(2*4*(1-1./q.^2)./(3*(2-2./q.^(4/3)./ft1))+4*(1-1./q.^4).*q.^2./(3*(2*q.^(2/3).*(2-1./q.^2)./ft1-2))),...
%     1./q,3*ft2./(2*4*(1-q.^2)./(3*(2-2*q.^(4/3)./ft2))+4*(1-q.^4)./q.^2./(3*(2./q.^(2/3).*(2-q.^2)./ft2-2)))); 
% grid