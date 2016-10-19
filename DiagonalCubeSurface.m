% program generates movie of the fact that a curved surface can be composed of straight lines

lam = 0:0.01:1;

[x,y] = meshgrid(lam,lam);

z = x.*(1-y) + (1-x).*y;

surf(x,y,z)
colormap copper
alpha(0.5)
shading interp
hold on
for j=1:11
    plot3(x(:,10*(j-1)+1),y(:,10*(j-1)+1),z(:,10*(j-1)+1),'r');
end
plot3([0 1],[0 0],[0 0],'b')
plot3([0 0],[0 1],[0 0],'b')
plot3([1 1],[0 1],[0 0],'b')
plot3([0 1],[1 1],[0 0],'b')
plot3([0 1],[0 0],[1 1],'b')
plot3([0 0],[0 1],[1 1],'b')
plot3([1 1],[0 1],[1 1],'b')
plot3([0 1],[1 1],[1 1],'b')
plot3([0 0],[0 0],[0 1],'b')
plot3([1 1],[0 0],[0 1],'b')
plot3([0 0],[1 1],[0 1],'b')
plot3([1 1],[1 1],[0 1],'b')
hold off
axis image
axis off
cameratoolbar
camlight

for j=1:60
    view(-37.5+(j-1)*6,30)
    eval(['print -dpng -r300 tmp' mint2str(j,2)])
end
