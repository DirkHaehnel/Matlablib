function Rotor3D(z, k, flag)

% Rotor3D(z, k) plots a 3D-plot of the spherical function z; k determines
% the number of rotations

if nargin<2 | isempty(k)
    k=1;
end

[n,m] = size(z);

theta = (0:pi/(n-1):pi)';
phi = 0:2*pi/(m-1):2*pi; 

figure(gcf)
if nargin<3 
    surf((sin(theta)*cos(phi)).*z,(sin(theta)*sin(phi)).*z,(cos(theta)*ones(size(phi))).*z,z)
else
    surf((sin(theta)*cos(phi)),(sin(theta)*sin(phi)),(cos(theta)*ones(size(phi))),z)
end
view([-37 10])
axis image
axis vis3d
shading interp
% axis off
whitebg([0.8 0.8 1])
tmp = axis;
tmp = tmp + max(abs(1.1*tmp-tmp))*sign(tmp);
hold on
plot3(tmp(1:2),[0 0],[0 0]);
plot3([0 0],tmp(3:4),[0 0]);
plot3([0 0],[0 0],tmp(5:6));
hold off
% colorbar
for j=1:round(m*k)
    camorbit(-360/(m-1),0)
    drawnow
end