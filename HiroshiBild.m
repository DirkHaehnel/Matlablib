set(0,'defaultLineLineWidth',1)

phi0 = pi/10; phi1 = pi/2;

theta=4/10*pi; phi=phi0;

pfeil([sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)]/2,[sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)]);
hold on
pfeil([sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)]/2,[0 0 0]);
psi = 0:pi/100:pi/2;
plot3(sin(psi)*cos(phi),sin(psi)*sin(phi),cos(psi))
psi = 0:pi/100:theta;
plot3(sin(psi)*cos(phi),sin(psi)*sin(phi),cos(psi),'b','linewidth',1.5)
plot3([0 cos(phi)],[0 sin(phi)],[0 0])

theta=1/6*pi; phi=phi1; 
pfeil([sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)]/2,[sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)]);
pfeil([sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)]/2,[0 0 0]);
psi = 0:pi/100:pi/2;
plot3(sin(psi)*cos(phi),sin(psi)*sin(phi),cos(psi))
psi = 0:pi/100:theta;
plot3(sin(psi)*cos(phi),sin(psi)*sin(phi),cos(psi),'b','linewidth',1.5)
plot3([0 cos(phi)],[0 sin(phi)],[0 0])

psi = 0:pi/100:2*pi;
plot3(cos(psi),sin(psi),psi*0)
psi = phi0:pi/100:phi1;
plot3(cos(psi),sin(psi),psi*0,'b','linewidth',1.5)
plot3([0 0],[0 0],[0 1.2])

hold off
axis image
