    % figure for 2color-2focus FCS

flag = 0;

close all
NA = 1.2;
chimax = asin(NA/1.333);
chi = (0:1e2)/100*chimax;
phi = pi/2 + (0:100)/200*2*pi;
foc = 200;
len = 2.1*foc*sin(chimax);

red = 0.8;
lighttransparency = 0.5;

d = 2;
len = 15*foc*sin(chimax);
height = 11*foc;
surf([0 len], [-len len], -foc*ones(2,2), zeros(2,2),'facealpha', 0.1) 
hold on
surf([0 0; len len], -len*ones(2,2), -foc*[1 d; 1 d], zeros(2,2),'facealpha', 0.3) 
surf(zeros(2,2), [-len -len; len len], -foc*[1 d; 1 d], zeros(2,2),'facealpha', 0.2)
theta0 = 0.8*chimax;
rad = foc*tan(theta0);
theta = asin(1.333*sin(theta0)/1.51);

% Light cone
surf([0; -rad]*cos(phi),[0; -rad]*sin(phi),[0; -foc]*ones(1,length(phi)),0.6*ones(2,length(phi)), 'facealpha', lighttransparency)
%surf([0; -3*rad]*cos(phi),[0; -3*rad]*sin(phi),[0; 3*foc]*ones(1,length(phi)),0.6*ones(2,length(phi)), 'facealpha', lighttransparency)

surf([-rad; -rad-tan(theta)*(d-1)*foc]*cos(phi),[-rad; -rad-tan(theta)*(d-1)*foc]*sin(phi),...
    [-foc; -d*foc]*ones(1,length(phi)),0.6*ones(2,length(phi)), 'facealpha', lighttransparency)
rad = rad+tan(theta)*(d-1)*foc;
surf([-rad; -rad-tan(theta0)*foc]*cos(phi), [-rad; -rad-tan(theta0)*foc]*sin(phi),...
    [-d*foc; -(d+1)*foc]*ones(1,length(phi)), 0.6*ones(2,length(phi)), 'facealpha', lighttransparency)
rad = rad+tan(theta0)*foc;
surf(-rad*[1;1]*cos(phi), -rad*[1;1]*sin(phi), [-(d+1)*foc*ones(1,length(phi)); -height-rad*sin(phi)], 0.6*ones(2,length(phi)), 'facealpha', lighttransparency)
surf(-rad*[1;1]*cos(phi), -rad*[1;1]*sin(phi), [-height-rad*sin(phi); -2*height*ones(1,length(phi))], 0.8*ones(2,length(phi)), 'facealpha', lighttransparency)
surf(-rad*[1;0]*cos(phi), -rad*[1;0]*sin(phi), [-2*height*ones(1,length(phi)); -3*height*ones(1,length(phi))], 0.8*ones(2,length(phi)), 'facealpha', lighttransparency)
surf(-rad*[0;1]*cos(phi), -rad*[0;1]*sin(phi), [-3*height*ones(1,length(phi)); -4*height*ones(1,length(phi))], 0.8*ones(2,length(phi)), 'facealpha', lighttransparency)
surf(-rad*[1;1]*cos(phi), -rad*[1;1]*sin(phi), [-4*height*ones(1,length(phi)); -5*height*ones(1,length(phi))], 0.8*ones(2,length(phi)), 'facealpha', lighttransparency)

len = 1.5*rad;
% Dichroic
surf(-rad*[1;1]*cos(phi), [-sin(phi)*rad; -height*ones(1,length(phi))], -height-rad*[1;1]*sin(phi), 0.6*ones(2,length(phi)), 'facealpha', lighttransparency)
surf([0 len], [-len len], -height-[1 1; -1 -1]*len, zeros(2,2),'facealpha', 0.1) 
surf([0 len], [-len len], -height-foc-[1 1; -1 -1]*len, zeros(2,2),'facealpha', 0.1)
surf([0 len; 0 len], -len*ones(2), -height-len+[-foc -foc; 0 0], zeros(2,2),'facealpha', 0.3) 
surf(zeros(2), -len*[1 -1; 1 -1], -height+[-foc-len -foc+len; -len +len], zeros(2,2),'facealpha', 0.3) 

% Laser beam combiner
surf([0 len], [-len len]-5*height, -height-[1 1; -1 -1]*len, zeros(2,2),'facealpha', 0.1) 
surf([0 len], [-len len]-5*height, -height-foc-[1 1; -1 -1]*len, zeros(2,2),'facealpha', 0.1)
surf([0 len; 0 len], -len*ones(2)-5*height, -height-len+[-foc -foc; 0 0], zeros(2,2),'facealpha', 0.3) 
surf(zeros(2), -len*[1 -1; 1 -1]-5*height, -height+[-foc-len -foc+len; -len +len], zeros(2,2),'facealpha', 0.3) 
surf(-rad*[1;1]*cos(phi), -5*height-rad*[1;1]*sin(phi), [-height-rad*sin(phi); -2.5*height*ones(1,length(phi))], 0.2*ones(2,length(phi)), 'facealpha', lighttransparency)
surf(-rad*[1;1]*cos(phi), [-5*height-rad*sin(phi); -6.5*height*ones(1,length(phi))], -height-rad*[1;1]*sin(phi), 0.5*ones(2,length(phi)), 'facealpha', lighttransparency)
surf(-rad*[1;1]*cos(phi), [-3*height*ones(1,length(phi)); -5*height-rad*sin(phi)], -height-rad*[1;1]*sin(phi), 0.6*ones(2,length(phi)), 'facealpha', lighttransparency)

% Phase modulator
%surf([0 len], [0 0]-3.5*height, -2*height-[1 1; -1 -1]*len, zeros(2,2),'facealpha', 0.1)
%surf([0:len/10:len]', [0 0]-3.5*height, -2*height-[1 1]'*linspace(1,-1,11)*len, zeros(2,11),'facealpha', 0.1)
surf([0 len], [0 0]-4*height, -height-[1 1; -1 -1]*len, zeros(2,2),'facealpha', 0.1)
surf([0:len/10:len]', [0 0]-4*height, -height-[1 1]'*linspace(1,-1,11)*len, zeros(2,11),'facealpha', 0.1)
surf(zeros(2,2), [-3.5 -3.5; -4 -4]*height, -height-[1 -1; 1 -1]*len, zeros(2,2),'facealpha', 0.2)
surf([0 0; len len], [-3.5 -4; -3.5 -4]*height, -height-[-1 -1; -1 -1]*len, zeros(2,2),'facealpha', 0.3)


% 1. Collimator Lens
rho = (0:50)'/50*rad;
lens = 1.5e3;
surf(-rho*cos(phi), -3*height-(sqrt(lens^2-rad^2)-sqrt(lens^2-rho.^2))*ones(1,length(phi)), -height-rho*sin(phi),  0*ones(length(rho),length(phi)), 'facealpha', 0.4)
surf(-rho*cos(phi), -3*height+(sqrt(lens^2-rad^2)-sqrt(lens^2-rho.^2))*ones(1,length(phi)), -height-rho*sin(phi), 0*ones(length(rho),length(phi)), 'facealpha', 0.4)
surf(-rad*[1;0]*cos(phi), -3*height+[0;5*foc]*ones(1,length(phi)), -height-rad*[1;0]*sin(phi), 0.5*ones(2,length(phi)), 'facealpha', lighttransparency)

% 2. Collimator Lens
rho = (0:50)'/50*rad;
lens = 1.5e3;
surf(-rho*cos(phi), -height-(sqrt(lens^2-rad^2)-sqrt(lens^2-rho.^2))*ones(1,length(phi)), -height-rho*sin(phi),  0*ones(length(rho),length(phi)), 'facealpha', 0.4)
surf(-rho*cos(phi), -height+(sqrt(lens^2-rad^2)-sqrt(lens^2-rho.^2))*ones(1,length(phi)), -height-rho*sin(phi), 0*ones(length(rho),length(phi)), 'facealpha', 0.4)
surf(-rad*[1;0]*cos(phi), -height-[0;5*foc]*ones(1,length(phi)), -height-rad*[1;0]*sin(phi), 0.5*ones(2,length(phi)), 'facealpha', lighttransparency)

% Waveguide
psi = (0:100)/100;
line(2.5*rad*sin(pi*psi).^4.*sin(2*pi*psi),-height-5*foc-(2*height-10*foc)*psi,-height-2.5*rad*sin(pi*psi).^4.*cos(2*pi*psi),'color','g')

% Objective
rho = (0:50)'/50*rad;
lens = 1.5e3;
surf(-rho*cos(phi), -rho*sin(phi), -(d+1)*foc-(sqrt(lens^2-rad^2)-sqrt(lens^2-rho.^2))*ones(1,length(phi)), 0*ones(length(rho),length(phi)), 'facealpha', 0.4)
surf(-rho*cos(phi), -rho*sin(phi), -(d+1)*foc+(sqrt(lens^2-rad^2)-sqrt(lens^2-rho.^2))*ones(1,length(phi)), 0*ones(length(rho),length(phi)), 'facealpha', 0.4)

% 1. Lens
rho = (0:50)'/50*rad;
lens = 1.5e3;
surf(-rho*cos(phi), -rho*sin(phi), -2*height-(sqrt(lens^2-rad^2)-sqrt(lens^2-rho.^2))*ones(1,length(phi)), 0*ones(length(rho),length(phi)), 'facealpha', 0.4)
surf(-rho*cos(phi), -rho*sin(phi), -2*height+(sqrt(lens^2-rad^2)-sqrt(lens^2-rho.^2))*ones(1,length(phi)), 0*ones(length(rho),length(phi)), 'facealpha', 0.4)

% 2. Lens
rho = (0:50)'/50*rad;
lens = 1.5e3;
surf(-rho*cos(phi), -rho*sin(phi), -4*height-(sqrt(lens^2-rad^2)-sqrt(lens^2-rho.^2))*ones(1,length(phi)), 0*ones(length(rho),length(phi)), 'facealpha', 0.4)
surf(-rho*cos(phi), -rho*sin(phi), -4*height+(sqrt(lens^2-rad^2)-sqrt(lens^2-rho.^2))*ones(1,length(phi)), 0*ones(length(rho),length(phi)), 'facealpha', 0.4)

% 3. Lens
rho = (0:50)'/50*rad;
lens = 1.5e3;
surf(-rho*cos(phi), -rho*sin(phi), -5*height-(sqrt(lens^2-rad^2)-sqrt(lens^2-rho.^2))*ones(1,length(phi)), 0*ones(length(rho),length(phi)), 'facealpha', 0.4)
surf(-rho*cos(phi), -rho*sin(phi), -5*height+(sqrt(lens^2-rad^2)-sqrt(lens^2-rho.^2))*ones(1,length(phi)), 0*ones(length(rho),length(phi)), 'facealpha', 0.4)

% Detectors
surf(-rad*[1;0]*cos(phi), -rad*[1;0]*sin(phi), -5*height-[0;5*foc]*ones(1,length(phi)), 0.8*ones(2,length(phi)), 'facealpha', lighttransparency)
%surf(-rad*[1;0]*cos(phi), height+[0;5*foc]*ones(1,length(phi)), -5*height-rad*[1;0]*sin(phi), 0.8*ones(2,length(phi)), 'facealpha', lighttransparency)
surf(-rad*[1;0]*cos(phi), -rad*[1;0]*sin(phi), -5*height-5*foc*ones(2,length(phi)), 0.*ones(2,length(phi)))
%surf(-rad*[0.4;0]*cos(phi), height+5*foc*ones(2,length(phi)), -6*height-rad*[0.4;0]*sin(phi), 0.*ones(2,length(phi)))

% Pinhole
surf(-[rad;100]*cos(phi), -[rad;100]*sin(phi), -foc/4-3*height*ones(2,length(phi)), 0.*ones(2,length(phi)), 'facealpha', 0.6)
surf(-[rad;100]*cos(phi), -[rad;100]*sin(phi), foc/4-3*height*ones(2,length(phi)), 0.*ones(2,length(phi)), 'facealpha', 0.6)
surf(-[rad;rad]*cos(phi), -[rad;rad]*sin(phi), [foc/4;-foc/4]*ones(1,length(phi))-3*height, 0.*ones(2,length(phi)), 'facealpha', 0.8)

% [x,y,z] = sphere(50);
% surf(30*x,30*y,30*z,zeros(size(x)));

hold off

caxis([0 1])
shading interp
axis image
axis off
view([-58 11])
material shiny
camlight
camlight