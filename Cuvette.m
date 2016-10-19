% drawing a cuvette measurement

phi = 0:pi/100:2*pi;
y = [-10 10];

a = 2;
h = 3;
d = 0.2;
for j=-1:2:1
    patch(ones(1,5)*a*j, [-a -a a a -a], [h -h -h h h], ones(1,5)*0.2, 'facealpha', 0.4-j/10);
    patch([-a -a a a -a], ones(1,5)*a*j, [h -h -h h h], ones(1,5)*0.2, 'facealpha', 0.4+j/10);    
end
patch([-a -a a a -a], [a -a -a a a], -h*ones(1,5), ones(1,5)*0.2, 'facealpha', 0.6);            
hold on
surf(cos(phi')*[1 1],ones(length(phi),1)*y,sin(phi')*[1 1],ones(length(phi),2)*0.6, 'facealpha', 0.7);
x = [0 10];
surf(ones(length(phi),1)*x, cos(phi')*[1 1], sin(phi')*[1 1],ones(length(phi),2)*1, 'facealpha', 0.6);
surf(cos(phi')*[0:0.1:1],ones(length(phi),11)*y(1),sin(phi')*[0:0.1:1],ones(length(phi),11)*0.6, 'facealpha', 0.7);




shading interp
axis image
axis off
camlight
view([60 10])
hold off
