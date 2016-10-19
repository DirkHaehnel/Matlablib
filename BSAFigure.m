figure%close

theta = (0:pi/300:pi*0.8)';
theta1 = (pi*0.8:pi/300:pi*0.9)';
phi = (0:300)/300*2*pi;
rad = 1;
hub = 0.;
thick = 0.1;

r = [rad*cos(theta),rad*sin(theta),hub*(theta/2-max(theta))/max(theta)];
tt = [-rad*sin(theta),rad*cos(theta),hub/2/max(theta)*ones(size(theta))];
tt = tt./(sqrt(sum(tt.^2,2))*[1 1 1]);
e1 = thick*[cos(theta),sin(theta),zeros(size(theta))];
e2 = cross(e1,tt,2);
e2 = thick*e2./(sqrt(sum(e2.^2,2))*[1 1 1]);

surf([r(1,1)*ones(size(phi)); r(1,1)*ones(size(phi)) + e1(1,1)*cos(phi) + e2(1,1)*sin(phi)],...
    [r(1,2)*ones(size(phi)); r(1,2)*ones(size(phi)) + e1(1,2)*cos(phi) + e2(1,2)*sin(phi)], ...
    [r(1,3)*ones(size(phi)); r(1,3)*ones(size(phi)) + e1(1,3)*cos(phi) + e2(1,3)*sin(phi)],...
    zeros(2,length(phi)));
hold on
surf(r(:,1)*ones(size(phi)) + e1(:,1)*cos(phi) + e2(:,1)*sin(phi), r(:,2)*ones(size(phi)) + e1(:,2)*cos(phi) + e2(:,2)*sin(phi), ...
    r(:,3)*ones(size(phi)) + e1(:,3)*cos(phi) + e2(:,3)*sin(phi), 0.8*theta/max(theta)*ones(size(phi)));
surf([r(end,1)*ones(size(phi)) + e1(end,1)*cos(phi) + e2(end,1)*sin(phi); r(end,1)*ones(size(phi)) + 1.2*(e1(end,1)*cos(phi) + e2(end,1)*sin(phi))],...
    [r(end,2)*ones(size(phi)) + e1(end,2)*cos(phi) + e2(end,2)*sin(phi); r(end,2)*ones(size(phi)) + 1.2*(e1(end,2)*cos(phi) + e2(end,2)*sin(phi))], ...
    [r(end,3)*ones(size(phi)) + e1(end,3)*cos(phi) + e2(end,3)*sin(phi); r(end,3)*ones(size(phi)) + 1.2*(e1(end,3)*cos(phi) + e2(end,3)*sin(phi))],...
    0.8*ones(2,length(phi)));
r = [rad*cos(theta1),rad*sin(theta1),hub*(theta1/2-max(theta))/max(theta)];
tt = [-rad*sin(theta1),rad*cos(theta1),hub/2/max(theta)*ones(size(theta1))];
tt = tt./(sqrt(sum(tt.^2,2))*[1 1 1]);
e1 = 1.2*thick*[cos(theta1),sin(theta1),zeros(size(theta1))];
e2 = cross(e1,tt,2);
e2 = 1.2*thick*e2./(sqrt(sum(e2.^2,2))*[1 1 1]);
diam = (max(theta1)-theta1)/(max(theta1)-min(theta1));
e1 = (diam*[1 1 1]).*e1;
e2 = (diam*[1 1 1]).*e2;
surf(r(:,1)*ones(size(phi)) + e1(:,1)*cos(phi) + e2(:,1)*sin(phi), r(:,2)*ones(size(phi)) + e1(:,2)*cos(phi) + e2(:,2)*sin(phi), ...
    r(:,3)*ones(size(phi)) + e1(:,3)*cos(phi) + e2(:,3)*sin(phi), 0.8*ones(length(theta1),length(phi)));
hold off
axis image
shading interp
colormap hot
caxis([-0.1 1])
%view([70-180 20])
view([70 20])
camlight(70-180, 20)
axis off
alpha(0.5)