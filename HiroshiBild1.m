set(0,'defaultLineLineWidth',1.5)

phi0 = 0; phi1 = pi/3; phi2 = 4*pi/5;
len = 5;
theta = 0; 
pfeil([0 0 0],len*[cos(theta)*cos(phi0),cos(theta)*sin(phi0),sin(theta)],[],[],[],0.7);
hold on
psi = phi2:pi/100:2*pi+phi1;
plot3(len*cos(psi),len*sin(psi),0*psi,'r')
psi = 0:pi/100:2*pi;
surf([0 len]'*cos(psi),[0 len]'*sin(psi),[0 0]'*psi,0.4*ones(2,length(psi),3),'facealpha',0.5)
pfeil([0 0 0],len*[cos(theta)*cos(phi1),cos(theta)*sin(phi1),sin(theta)],[],[],[],0.7);
theta = pi/5;
pfeil([0 0 0],len*[cos(theta)*cos(phi2),cos(theta)*sin(phi2),sin(theta)],[],[],[],0.7);
plot3(cos(phi2)*[0 len],sin(phi2)*[0 len],[0 0],'r')
psi = 0:pi/300:theta;
plot3(len*cos(phi2)*cos(psi),len*sin(phi2)*cos(psi),len*sin(psi),'b','linewidth',2.5)
psi = phi1:pi/300:phi2;
plot3(len*cos(psi),len*sin(psi),0*psi,'b','linewidth',2.5)



hold off
axis image
