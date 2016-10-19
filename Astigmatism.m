% program Astigmatism for the calculation of astigmatic laser focus

close all
clear all
zfield = [0.025 4];
NA = 1.2;
delta = 0.04;
n0 = 1.333;
n = 1.333;
n1 = 1.333;
d0 = [];
d = 4;
d1 = [];
lamex = 0.67;
ww = 3.5e3;
wd = 3e3;
focpos = 2; 
zetav = 0:0.1:0.4;

load M1
[tx,ty] = meshgrid(0:size(res,2)-1,0:size(res,1)-1);
tx=(tx-max(tx(:))/2)*delta; ty=(ty-max(ty(:))/2)*delta;
r = sqrt(tx.^2+ty.^2);
p = angle(tx+i*ty)-25/180*pi;
rhofield = [-0.025 1.2*max(r(:))];

phi = 0:pi/50:2*pi;
row = ones(size(phi));
cnt = 1;
for jzeta=1:length(zetav)
    zeta = zetav(jzeta)*2*pi*ww^2/NA^2/wd^2; 
    w0 = ww/sqrt(1+zeta^2);
    z0 = pi*w0^2*zeta/lamex;
    over = [NA*3e3  -z0 z0  w0 w0];
    [fx0, fx2, fz, rho, z, exc, excx, excz] = NewExc(rhofield, zfield, NA, n0, n, n1, d0, d, d1, lamex, over, focpos, [1.52 0], lamex/0.05);
    ind = 1:size(z,2);
    for j = 1:3
        s = ind(abs(z(1,:)-focpos+0.-0.25*(j-2))<eps);
        if size(fx0,3)==1
            int(:,:,j) = abs(interp1(rho(:,1),fx0(:,s),r,'cubic') + interp1(rho(:,1),fx2(:,s),r,'cubic').*cos(2*p)).^2;
            int(:,:,j) = int(:,:,j) + abs(interp1(rho(:,1),fx2(:,s),r,'cubic').*sin(2*p)).^2;
        else
            intx(:,:,j) = interp1(rho(:,1),fx0(:,s,1),r,'cubic'); 
            inty(:,:,j) = 0*intx(:,:,j);
            intz(:,:,j) = inty(:,:,j);
            for k=2:4
                intx(:,:,j) = intx(:,:,j) + interp1(rho(:,1),fx0(:,s,k),r,'cubic').*cos(2*(k-1)*p);
                inty(:,:,j) = inty(:,:,j) + interp1(rho(:,1),fx2(:,s,k),r,'cubic').*sin(2*(k-1)*p);
                intz(:,:,j) = intz(:,:,j) + interp1(rho(:,1),fz(:,s,k),r,'cubic').*cos(2*(k-3)*p);
            end
            int(:,:,j) = 1^2*abs(intx(:,:,j)).^2 + (1-1^2)*abs(inty(:,:,j)).^2; % + abs(intz(:,:,j)).^2;
        end
    end
    for j=1:3
        s = ind(z(1,:)==focpos+0.25*(j-2));        
        subplot(5,length(zetav)+1,(j-1)*(length(zetav)+1)+jzeta)
        pcolor(tx,ty,int(:,:,j)/max(int(:)));
        caxis([0 1]);
        set(gca,'fontsize',10);
        if j==1
            title([num2str(zetav(jzeta)) '\times2\pi'],'fontsize',12)
        end
        if jzeta==1
            text(-1.2*max(tx(:)),0,[num2str(1e3*(z(1,s)-focpos)) ' nm'],'fontsize',12,'horizontalalignment','right')
        end
        axis image; shading interp
        axis off; colormap hot
        drawnow
        cnt = cnt+1;
    end
end
for j=1:3
    subplot(5,length(zetav)+1,j*length(zetav)+j)
    pcolor(tx,ty,res(:,:,4-j)); %/max(max(max(res(:,:,1:3)))));
    %caxis([0 1]);        
    axis image; shading interp
    axis off; colormap hot
    if j==1 title('exp.','fontsize',12); end
    drawnow
end    

