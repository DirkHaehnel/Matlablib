% figure for astigmatism paper

close all
NA = 1.2;
chimax = asin(NA/1.333);
chi = (0:1e2)/100*chimax;
phi = pi/2 + (0:100)/200*2*pi;
foc = 200;
len = 2.1*foc*sin(chimax);

d = 1;
dd = 25;
plot3([0 0], [0 0], [-d*foc 0.5*foc], 'linewidth', 1, 'color', [0,0,0])
hold on
theta0 = 0.8*chimax;
theta = asin(1.333*sin(theta0)/1.51);
% surf([zeros(1,length(phi)); tan(theta0)*(-foc+dd*cos(2*phi)).*cos(phi)],...
%     [zeros(1,length(phi)); tan(theta0)*(-foc+dd*cos(2*phi)).*sin(phi)],...
%     [dd*cos(2*phi); -foc*ones(1,length(phi))],0.6*ones(2,length(phi)), 'facealpha', 0.4)
[x,y,z] = sphere(50);
surf(3*x,3*y,3*z+dd,zeros(size(x)));
surf(3*x,3*y,3*z-dd,zeros(size(x)));        
[x,y,z] = cylinder;
surf(3*x,3*y,2*(z-0.5)*dd,zeros(size(x)));                
len = 120;
surf(-len*sin(chi')*cos(phi),len*sin(chi')*sin(phi),-len*cos(chi')*ones(size(phi)),0.*ones(length(chi), length(phi)), 'facealpha', 0.1)
surf(-len*sin(chi')*cos(phi),len*sin(chi')*sin(phi),-sin(chi').^2*dd*cos(2*phi)-len*cos(chi')*ones(size(phi)),0.*ones(length(chi), length(phi)), 'facealpha', 0.3)
chi = max(chi);
% surf([zeros(1,length(phi)); -len*sin(chi')*cos(phi)],...
%     [zeros(1,length(phi)); len*sin(chi')*sin(phi)],...
%     [dd*cos(2*phi); -sin(chi').^2*dd*cos(2*phi)-len*cos(chi')*ones(size(phi))],0.6*ones(2,length(phi)), 'facealpha', 0.4)

phi = pi/2;
plot3([0; -len*sin(chi)*cos(phi)],[0; len*sin(chi)*sin(phi)], [dd*cos(2*phi); -sin(chi).^2*dd*cos(2*phi)-len*cos(chi')*ones(size(phi))])
phi = pi;
plot3([0; -len*sin(chi)*cos(phi)],[0; len*sin(chi)*sin(phi)], [dd*cos(2*phi); -sin(chi).^2*dd*cos(2*phi)-len*cos(chi')*ones(size(phi))])
phi = 3*pi/2;
plot3([0; -len*sin(chi)*cos(phi)],[0; len*sin(chi)*sin(phi)], [dd*cos(2*phi); -sin(chi).^2*dd*cos(2*phi)-len*cos(chi')*ones(size(phi))])


hold off



caxis([0 1])
shading interp
axis image
axis off
view([-58 11])
