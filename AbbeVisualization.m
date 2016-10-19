close all

if 0 % lateral resolution
    my = 50;
    mx = 250;
    al=30/180*pi;

    [xx,yy] = meshgrid(-mx:0.5:0,-mx*tan(al):0.5:my+mx*tan(al));

    mask1 = yy<=my+tan(al)*xx & yy>=tan(al)*xx;
    mask2 = yy<=my-tan(al)*xx & yy>=-tan(al)*xx;

    pcolor(xx,yy,(mask1 | mask2).*(2+real(mask1.*exp(i*(cos(al)*xx+sin(al)*yy))+mask2.*exp(i*(cos(al)*xx-sin(al)*yy)))));
    axis equal; colormap(flipud(hot)); shading interp
    axis off

    % axis([-50 0 -40 my+40])
end

if 0 % axial resolution
    my = 50;
    mx = 250;
    al=30/180*pi;

    [xx,yy] = meshgrid(-mx:0.5:mx/2,-mx*tan(al):0.5:my+mx*tan(al));

    mask1 = yy<=my+tan(al)*xx & yy>=tan(al)*xx;
    mask2 = yy<=(my+my*cos(al))/2 & yy>=(my-my*cos(al))/2;

    pcolor(xx,yy,(mask1 | mask2).*(2+real(mask1.*exp(i*(cos(al)*xx+sin(al)*yy))+mask2.*exp(i*xx))));
    axis equal; colormap(flipud(hot)); shading interp
    axis off

    % axis([-50 0 -40 my+40])
end

if 0 % focus generation
    my = 50;
    mx = 50;
    al = (0:0.1:60)/180*pi;
    cnt = 1;

    [xx,yy] = meshgrid(-mx:0.5:mx/2,-my/2-mx*tan(max(al)):0.5:my/2+mx*tan(max(al)));
    feld = 0*xx; focus = 0*xx;
    mask0 = xx>-25;
    for j=1:length(al)
        mask = xx.*sin(al(j))+yy.*cos(al(j)) < 10 & xx.*sin(al(j))-yy.*cos(al(j)) < 10;
        focus = focus + exp(i*(cos(al(j))*xx-sin(al(j))*yy)) + exp(i*(cos(al(j))*xx+sin(al(j))*yy));
        feld = focus;
        feld(~mask & ~mask0) = nan;
        feld(mask0) = 0;
        if rem(j-1,10)==0
            tmp = abs(focus).^2/max(abs(focus(:)).^2)*max(real(focus(:)))*2 + min(real(focus(:)));
            tmp(~mask0) = nan;
            pcolor([xx 3/2*mx+xx],[yy yy],[real(feld) + mask0.*real(focus) tmp]);
            axis equal; colormap(hot); shading interp
            axis off
            pause(0.01); drawnow
            eval(['print -dpng -r300 AbbeFocus' mint2str(cnt,2)]); cnt = cnt+1;
        end
    end
end

if 0 % full laser focus
    rhofield = [0 0.75];
    zfield = [-1 1];
    NA = 1.2;
    fd = 3e3;
    n0 = 1.33;
    n = 1.33;
    n1 = 1.33;
    d0 = [];
    d = 0;
    d1 = [];
    lamex = 0.514;
    over = 5e3;
    focpos = 0;
    atf = [];
    resolution = [25 25];
    [fxc, fxs, fyc, fys, fzc, fzs, rho, z] = GaussExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, atf, resolution);
    FocusImage2D(rho,z,cat(4,cat(3,fxc,fxs),cat(3,fyc,fys),cat(3,fzc,fzs)),0,'horizontal')
    
    figure
    fxc=fxc+fxc(:,end:-1:1,:);fxs=fxs+fxs(:,end:-1:1,:);fyc=fyc+fyc(:,end:-1:1,:);fys=fys+fys(:,end:-1:1,:);fzc=fzc+fzc(:,end:-1:1,:);fzs=fzs+fzs(:,end:-1:1,:);
    FocusImage2D(rho,z,cat(4,cat(3,fxc,fxs),cat(3,fyc,fys),cat(3,fzc,fzs)),0,'horizontal')
end

if 0 % lateral resolution
    len = 200;
    [x,y] = meshgrid(-len:0.2:0,-len:0.2:len);
    sze = length(x(1,:));
    for j=3.5:-0.1:0.6
        z = 2*pi*j;
        ind = sqrt(x.^2+y.^2)>len;
        feld = real(exp(i*sqrt(x.^2+(y+z/2).^2))+exp(i*sqrt(x.^2+(y-z/2).^2)));
        feld(ind) = 2.1;
        mim(feld)
        phi = asin(pi/z);
        line([sze sze*(1-cos(phi))],[sze sze*(1+sin(phi))])
        line([sze sze*(1-cos(phi))],[sze sze*(1-sin(phi))])
        drawnow
        eval(['print -dpng -r300 AbbeInterference' mint2str(36-10*j,2)])
    end
end

if 0 % axial resolution (crap)
    len = 200;
    [x,y] = meshgrid(-len:0.2:0,-len:0.2:len);
    sze = length(x(1,:));
    theta = asin(1.2/1.33);
    z = pi/(1-cos(theta));
    ind = sqrt(x.^2+y.^2)>len;
    %feld = unwrap(angle(exp(i*sqrt((x-z/2).^2+y.^2))+exp(i*sqrt((x+z/2).^2+y.^2))));
    feld = unwrap(angle(1+exp(i*sqrt((x+z/2).^2+y.^2))./exp(i*sqrt((x-z/2).^2+y.^2))));
    feld(ind) = 2.1;
    mim(feld)
    phi = acos(1-2*pi/z);
    line([sze sze*(1-cos(theta))],[sze sze*(1+sin(theta))])
    line([sze sze*(1-cos(theta))],[sze sze*(1-sin(theta))])
    drawnow
    %eval(['print -dpng -r300 AbbeAxialInterference' mint2str(10*j,2)])
end

if 1 % Structured illumination
    [xx, yy] = meshgrid(0:0.1:6*pi,0:0.1:10);
    pcolor(xx,yy,abs(exp(1i*xx)+exp(-1i*xx)).^2)
    colormap(flipud(hot)); shading interp
    axis off
    axis image

    figure
    [xx, yy] = meshgrid(pi/4:0.1:6.25*pi,0:0.1:10);
    pcolor(xx,yy,abs(exp(1i*xx)+exp(-1i*xx)).^2)
    colormap(flipud(hot)); shading interp
    axis off
    axis image

    figure
    pcolor(xx,yy,0.5*ones(size(xx)))    
    colormap(flipud(hot)); shading interp
    axis off
    axis image
    
    
end
