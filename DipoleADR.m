tst = [];
lw = 1;

% alv = 45/180*pi; bev = 225/180*pi;  

%alv = [45 45]; bev = [0 360]; nn = 36;
%alv = [0 90]; bev = [0 0]; nn = 9;
%alv = [90 90]; bev = [0 360]; nn = 36;

r = ones(37,1)*[1 0 1]/2 + cos((0:10:360)'/180*pi)*[1 0 -1]/2 - sin((0:10:360)'/180*pi)*[0 -1 0]/sqrt(2);
alv = acos(r(:,3));
bev = pi + atan(r(:,2)./r(:,1));

% if ~(diff(alv)==0)
%     alv = (alv(1):diff(alv)/nn:alv(2))/180*pi;
% else
%     alv = alv(1)*ones(1,nn+1)/180*pi;
% end
% if ~(diff(bev)==0)
%     bev = pi + (bev(1):diff(bev)/nn:bev(2))/180*pi;
% else
%     bev = pi + bev(1)*ones(1,nn+1)/180*pi;
% end

for j=1:length(alv)

    al = alv(j);
    be = bev(j);

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

    phi=(0:200)'/100*pi;
    col = ones(size(col));
    [v,pc,ps] = DipoleL(psi,0,n2,n1,n1,[],0,[]);
    plot3((cos(phi)*sin(psi)).*(abs(cos(al)*col*v+sin(al)*cos(phi-be)*pc).^2+abs(sin(al)*sin(phi-be)*ps).^2),...
        (sin(phi)*sin(psi)).*(abs(cos(al)*col*v+sin(al)*cos(phi-be)*pc).^2+abs(sin(al)*sin(phi-be)*ps).^2),...
        -(col*cos(psi)).*(abs(cos(al)*col*v+sin(al)*cos(phi-be)*pc).^2+abs(sin(al)*sin(phi-be)*ps).^2),'b','linewidth',lw)

    colormap jet
    axis image
    shading flat
    alpha(0.5)

    [intx, inty, intz, rho, phi] = SEPDipole([0 2], 0, 1.24, n2, n1, n1, [], 0, [], 0.67, 100, 1, [], [], [al be]);
    rho = rho/50;
    surf(cos(phi)*rho,sin(phi)*rho,-10*ones(length(phi),length(rho)),intx/max(intx(:))); axis image; shading interp

    light 
    
    view([20 30]);
    tst = [tst; axis];
    axis([-9.5 9.5 -4 9.5 -11 2])
    
	hold off

    eval(['print -dpng -r300 DipoleADR' mint2str(j,2)])
%     if j<=nn
%         eval(['print -dpng -r300 DipoleADR' mint2str(2*nn+2-j,2)])
%     end        

end
