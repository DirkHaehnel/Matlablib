theta = (0:pi/100:pi)';
phi = 0:pi/100:2*pi;
a = permute(cat(3,zeros(3),0.15*eye(3)),[3,1,2]);

subplot(3,5,3)
f = real(SphericalHarmonic(0,0,theta,phi));
surf((sin(theta)*cos(phi)).*abs(f).^2,(sin(theta)*sin(phi)).*abs(f).^2,(cos(theta)*ones(size(phi))).*abs(f).^2,abs(f).^2);
alpha(0.5)
camlight
line(a(:,:,1),a(:,:,2),a(:,:,3),'color',[0.5 0.3 1])
view([137.5 30])
axis off
axis image

subplot(3,5,3)
f = real(SphericalHarmonic(0,0,theta,phi));
surf((sin(theta)*cos(phi)).*abs(f).^2,(sin(theta)*sin(phi)).*abs(f).^2,(cos(theta)*ones(size(phi))).*abs(f).^2,abs(f).^2);
line(a(:,:,1),a(:,:,2),a(:,:,3),'color',[0.5 0.3 1])
view([137.5 30])
axis off
axis image


f = real(SphericalHarmonic(1,0,theta,phi));
surf((sin(theta)*cos(phi)).*abs(f).^2,(sin(theta)*sin(phi)).*abs(f).^2,(cos(theta)*ones(size(phi))).*abs(f).^2,abs(f).^2);
a = permute(cat(3,zeros(3),0.15*eye(3)),[3,1,2]);
line(a(:,:,1),a(:,:,2),a(:,:,3),'color',[0.5 0.3 1])
view([137.5 30])
axis off
axis image

