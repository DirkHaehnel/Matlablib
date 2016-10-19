r = [1 2 3 5 6];
phi = (0:pi/100:3/2*pi)';
z = 0:0.1:max(r)*(3+sqrt(5))/2;

al = 1./(0.5+0.5*(1:length(r)));
co = [ones(length(r),1) (0:length(r)-1)'/length(r) zeros(length(r),1)];

j = 1;
surf(r(j)*cos(phi)*ones(size(z)),r(j)*sin(phi)*ones(size(z)),ones(size(phi))*z,'FaceColor',co(j,:),'EdgeColor','none','FaceAlpha',al(j));
hold on
for j=1:length(r)
    surf(r(j)*cos(phi)*ones(size(z)),r(j)*sin(phi)*ones(size(z)),ones(size(phi))*z,'FaceColor',co(j,:),'EdgeColor','none','FaceAlpha',al(j));
end
hold off
axis image
view([30 30])
axis off