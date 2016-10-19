tst = [];
lw = 1;

al = -45/180*pi;
be = 30/180*pi;;

theta = (0:300)/300*pi/2;
phi=(0:200)'/200*pi;
col = ones(size(phi));

n1 = 1.0;
n2 = 1.52;
psi = asin(n1/n2);
theta = sort([psi theta]);

[v,pc,ps] = DipoleL(theta,0,n2,n1,n1,[],0,[]);
v = v.';
pc = pc.';
ps = ps.';

strt = [0 0 0];
Pfeil(strt+[0 0 0],strt+1.5*[cos(be)*sin(al) sin(be)*sin(al) cos(al)])

hold on

line([strt(1) strt(1)-3],[strt(2) strt(2)],[strt(3) strt(3)],'color','k','linewidth',1)
line([strt(1) strt(1)],[strt(2) strt(2)-1.5],[strt(3) strt(3)],'color','k','linewidth',1)
line([strt(1) strt(1)],[strt(2) strt(2)],[strt(3) strt(3)+1.5],'color','k','linewidth',1)

text(1,0,2,['\theta = ' int2str(al/pi*180) '°, \phi = ' int2str(be/pi*180-180) '°'])

surf((cos(phi)*sin(theta)).*(abs(cos(al)*col*v+sin(al)*cos(phi-be)*pc).^2+abs(sin(al)*sin(phi-be)*ps).^2),...
    (sin(phi)*sin(theta)).*(abs(cos(al)*col*v+sin(al)*cos(phi-be)*pc).^2+abs(sin(al)*sin(phi-be)*ps).^2),...
    -(col*cos(theta)).*(abs(cos(al)*col*v+sin(al)*cos(phi-be)*pc).^2+abs(sin(al)*sin(phi-be)*ps).^2))

col = 1;
phi = 0;
plot3((cos(phi)*sin(theta)).*(abs(cos(al)*col*v+sin(al)*cos(phi-be)*pc).^2+abs(sin(al)*sin(phi-be)*ps).^2),...
    (sin(phi)*sin(theta)).*(abs(cos(al)*col*v+sin(al)*cos(phi-be)*pc).^2+abs(sin(al)*sin(phi-be)*ps).^2),...
    -(col*cos(theta)).*(abs(cos(al)*col*v+sin(al)*cos(phi-be)*pc).^2+abs(sin(al)*sin(phi-be)*ps).^2),'r','linewidth',lw)

phi = pi;
plot3((cos(phi)*sin(theta)).*(abs(cos(al)*col*v+sin(al)*cos(phi-be)*pc).^2+abs(sin(al)*sin(phi-be)*ps).^2),...
    (sin(phi)*sin(theta)).*(abs(cos(al)*col*v+sin(al)*cos(phi-be)*pc).^2+abs(sin(al)*sin(phi-be)*ps).^2),...
    -(col*cos(theta)).*(abs(cos(al)*col*v+sin(al)*cos(phi-be)*pc).^2+abs(sin(al)*sin(phi-be)*ps).^2),'r','linewidth',lw)

phi=(0:100)'/100*pi;
col = ones(size(col));
[v,pc,ps] = DipoleL(psi,0,n2,n1,n1,[],0,[]);
plot3((cos(phi)*sin(psi)).*(abs(cos(al)*col*v+sin(al)*cos(phi-be)*pc).^2+abs(sin(al)*sin(phi-be)*ps).^2),...
    (sin(phi)*sin(psi)).*(abs(cos(al)*col*v+sin(al)*cos(phi-be)*pc).^2+abs(sin(al)*sin(phi-be)*ps).^2),...
    -(col*cos(psi)).*(abs(cos(al)*col*v+sin(al)*cos(phi-be)*pc).^2+abs(sin(al)*sin(phi-be)*ps).^2),'b','linewidth',lw)

colormap jet
axis image
shading flat
alpha(0.5)

[intx, inty, intz, rho, phi] = SEPDipole([0 2], 0, 1.24, n2, n1, n1, [], 0, [], 0.67, 100, 1, [], [], [], [be al]);
rho = rho/50;
surf(cos(phi)*rho,sin(phi)*rho,-10*ones(length(phi),length(rho)),intx/max(intx(:))); axis image; shading interp

surf([-9.5 9.5 9.5 -9.5],[0 0 9.5 9.5],[0 0 0 0],[[0.2 0.2 0.2 0.2];zeros(2,4)]);%,'alpha',0.2); axis image; shading interp

light

view([20 30]);
axis([-9.5 9.5 -4 9.5 -11 2])

hold off

