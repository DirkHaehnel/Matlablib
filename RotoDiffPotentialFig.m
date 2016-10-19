% RotoDiffPotentialFig

[x,y,z]=ellipsoid(0,0,0,1,1,0.5,100);
surf(x,y,z,sqrt(x.^2+y.^2+z.^2))
shading interp
axis off
camlight
axis image
colomap copper
material metal
alpha(0.5)