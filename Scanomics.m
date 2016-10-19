close
z=(-3:0.1:3)';
phi=0:pi/20:2*pi;
num = 1e3;
[x, y] = meshgrid(-2.3:0.1:2.3,-2.3:0.1:2.3);

surf(abs(z)/1.4*cos(phi),abs(z)/1.4*sin(phi),z*ones(size(phi)),ones(size(z))*ones(size(phi)));
hold
c = 1.5*exp(-(x.^2+y.^2));
c(end-1:end,end-1:end) = 0;
surf(x,y,zeros(size(x)),0.5+exp(-(x.^2+y.^2)));
hold
axis image; axis off; view([-31,20]); colormap copper; shading interp; material shiny; lighting phong; light; light

