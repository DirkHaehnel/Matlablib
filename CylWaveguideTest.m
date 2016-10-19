% CylWaveguideTest

n1 = 1.52; n2 = 1.46; rad = 2*pi; m = 1;

[qp, coef] = CylWaveguideMode([n1 n2],rad,m);

rr = (0.5:1e3)/1e3*5*rad;

k = 1;
phi = (0:1e2)'/1e2*2*pi;
[fr, ff, fz] = CylField([n1 n2],rad,m,qp(k),coef(:,k),rr,phi);

% subplot(131); surf(cos(phi)*rr/2/pi,sin(phi)*rr/2/pi,abs(fr).^2); shading interp; view([-36 10]); grid off; light
% subplot(132); surf(cos(phi)*rr/2/pi,sin(phi)*rr/2/pi,abs(ff).^2); shading interp; view([-36 10]); grid off; light
% subplot(133); surf(cos(phi)*rr/2/pi,sin(phi)*rr/2/pi,abs(fz).^2); shading interp; view([-36 10]); grid off; light
% colormap hot
% 
subplot(131); pcolor(cos(phi)*rr/2/pi,sin(phi)*rr/2/pi,abs(fr).^2); axis image; shading interp; colorbar('h'); 
subplot(132); pcolor(cos(phi)*rr/2/pi,sin(phi)*rr/2/pi,abs(ff).^2); axis image; shading interp; colorbar('h'); 
subplot(133); pcolor(cos(phi)*rr/2/pi,sin(phi)*rr/2/pi,abs(fz).^2); axis image; shading interp; colorbar('h'); 
colormap hot

figure

plot(rr/2/pi,imag(fr(1,:))+real(fr(1,:)),rr/2/pi,imag(ff(1,:))+real(ff(1,:)),rr/2/pi,imag(fz(1,:))+real(fz(1,:)))
ax = axis;
line([1 1]*rad/2/pi,ax(3:4),'linestyle',':');

%quiver(cos(phi)*rr/2/pi,sin(phi)*rr/2/pi,imag(cr.*fr-sr.*ff),imag(sr.*fr+cr.*ff),0.05); axis image