% program for calculating defoc excitation scan images of a dipole emitter
% within a diamond film
%
% all length values are in units of micrometers!!

%%
lamem = 0.76; % excitation wavelength
resolution = 50;
molpos = 0.01:0.01:1; % molecule's distance in diamond layer from surface
mag = 100;
pixel = 6.45; % CCD pixel size
nn = 40;
rhofield = [-lamem/resolution(1)/2 1.5*nn*pixel/mag];
NA = 1.3; %0.95; % N.A. of objective
fd = 3e3; % foccal distance of objective (Olympus)
n = 2.417; % diamond ref. index
n0 = 1.52; %1.0; % immersion ref. index
n1 = n;
d0 = [];
d1 = [];
focpos = 1; % vector of defoc values (should be 48 values!!, currently in steps of 50 nm)
% out of plane rotation of dipole axis
theta = [90 0]/180*pi;
% in-plane rotation of dipole axis
phi = zeros(size(theta)); 
%%

%%
[x,y] = meshgrid(-nn:nn,-nn:nn);
p = angle(x+1i*y);
r = sqrt(x.^2+y.^2);
mask = zeros(2*nn+1,2*nn+1,length(molpos),length(theta));
for jz=1:length(molpos)
    [intx, inty, intz, rho, ~, fxx0, fxx2, fxz, byx0, byx2, byz] = ...
        SEPDipole(rhofield, molpos(jz), NA, n0, n, n1, d0, molpos(jz), d1, lamem, mag, focpos);
    for cnt=1:length(theta)
        al = theta(cnt);
        be = phi(cnt);
        mask(:,:,jz,cnt) = SEPImage(al,be,nn,pixel,rho,fxx0,fxx2,fxz,byx0,byx2,byz);
    end
end
for j=1:10 sx{j} = [int2str(molpos(j)*1e3) ' nm']; end
for j=1:10 sy{j} = ['+' int2str((molpos((j-1)*10+1)-molpos(1))*1e3) ' nm']; end
sy{1} = '';

% shows the defoc images for different defoc values for an inclinatin angle
% of theta(11); change '11' to other values to see result for other inclination angles
CombineImages(mask(:,:,:,1),10,10,'scale',sx,sy);
set(gca,'Position',[0.1 0.1 0.8 0.8])
bla = get(gca,'children');
for j=1:length(bla) set(bla(j),'fontsize',12); end
%%

%%
close
phi = 0:pi/5e3:pi/2;
[v,pc,ps] = DipoleL(phi,0.01/lamem*2*pi,n0,n,n1,d0,0.01/lamem*2*pi,d1);
[v1,pc1,ps1] = DipoleL(phi,0,n1,n,n0,d0,0.01/lamem*2*pi,d1);
plot([0 4],[0 4/tan(asin(n0/n))],'color',[0.8 0.8 0.8])
hold on
plot([0 -4],[0 4/tan(asin(n0/n))],'color',[0.8 0.8 0.8])
plot(sin(phi)'.*abs(v).^2,-cos(phi)'.*abs(v).^2,...
    -sin(phi)'.*abs(v).^2,-cos(phi)'.*abs(v).^2,'r',...
    sin(phi)'.*(abs(pc).^2+abs(ps).^2)/2,-cos(phi)'.*(abs(pc).^2+abs(ps).^2)/2,'b',...
    -sin(phi)'.*(abs(pc).^2+abs(ps).^2)/2,-cos(phi)'.*(abs(pc).^2+abs(ps).^2)/2,'b',...
    sin(phi)'.*abs(v1).^2,cos(phi)'.*abs(v1).^2,'r',...
    -sin(phi)'.*abs(v1).^2,cos(phi)'.*abs(v1).^2,'r',...
    sin(phi)'.*(abs(pc1).^2+abs(ps1).^2)/2,cos(phi)'.*(abs(pc1).^2+abs(ps1).^2)/2,'b',...
    -sin(phi)'.*(abs(pc1).^2+abs(ps1).^2)/2,cos(phi)'.*(abs(pc1).^2+abs(ps1).^2)/2,'b');
plot([-4 4],[0 0],'color',[0.8 0.8 0.8])
hold off
text(2,0.3,'diamond','FontSize',20)
text(2,-0.3,'oil (\itn\rm = 1.52)','FontSize',20)
%%

%%
close
molpos = 0.001:0.001:0.1;
for jz=1:length(molpos)
    [lvd,lvu,lpd,lpu,qvd,qvu,qpd,qpu,qv,qp] = LifetimeL(molpos(jz)/lamem*2*pi,n0,n,n1,d0,molpos(jz)/lamem*2*pi,d1);
    taup(jz) = 4/3*n/(qpd+qpu);
    tauv(jz) = 4/3*n/(qvd+qvu);
end
plot(molpos*1e3,[taup;tauv])
xlabel('distance from surface (nm)')
ylabel('rel. lifetime')
legend({'horizontal dipole','vertical dipole'})
%%

%%
lamem = 0.76; % excitation wavelength
resolution = 50;
molpos = 0.01:0.01:1; % molecule's distance in diamond layer from surface
mag = 100;
pixel = 6.45; % CCD pixel size
nn = 40;
rhofield = [-lamem/resolution(1)/2 1.5*nn*pixel/mag];
fd = 3e3; % foccal distance of objective (Olympus)
n = 2.417; % diamond ref. index
n0 = 1.52; %1.0; % immersion ref. index
n1 = n;
d0 = [];
d1 = [];
NA = 1.3/n0*n; 
focpos = 1; % vector of defoc values (should be 48 values!!, currently in steps of 50 nm)
% out of plane rotation of dipole axis
theta = [90 0]/180*pi;
% in-plane rotation of dipole axis
phi = zeros(size(theta)); 

[x,y] = meshgrid(-nn:nn,-nn:nn);
p = angle(x+1i*y);
r = sqrt(x.^2+y.^2);
mask = zeros(2*nn+1,2*nn+1,length(molpos),length(theta));
for jz=1:length(molpos)
    [intx, inty, intz, rho, ~, fxx0, fxx2, fxz, byx0, byx2, byz] = ...
        SEPDipole(rhofield, 0, NA, n1, n, n0, [], molpos(jz), [], lamem, mag, focpos);
    for cnt=1:length(theta)
        al = theta(cnt);
        be = phi(cnt);
        mask(:,:,jz,cnt) = SEPImage(al,be,nn,pixel,rho,fxx0,fxx2,fxz,byx0,byx2,byz);
    end
end
for j=1:10 sx{j} = [int2str(molpos(j)*1e3) ' nm']; end
for j=1:10 sy{j} = ['+' int2str((molpos((j-1)*10+1)-molpos(1))*1e3) ' nm']; end
sy{1} = '';

% shows the defoc images for different defoc values for an inclinatin angle
% of theta(11); change '11' to other values to see result for other inclination angles
CombineImages(mask(:,:,:,1),10,10,'scale',sx,sy);
set(gca,'Position',[0.1 0.1 0.8 0.8])
bla = get(gca,'children');
for j=1:length(bla) set(bla(j),'fontsize',12); end
%%