close all
pixel_size = 0.06;

if 0
    [fname, dname] = uigetfile;
else
    dname = 'C:\Data\16042012\stack';
    fname = dir([dname '*.mat']);
end


if isstruct(fname)
    for j=1:length(fname)
        load([dname fname(j).name])
        imraw = reshape(results.SHR,round([results.Info.xCCD results.Info.yCCD results.Info.xPixel results.Info.yPixel]));
        
        mim(sum(sum(imraw,3),4))
        eval(['print -dpng -r300 ''' dname fname(j).name(1:end-4) '_CCD'''])
        
        step_size = (results.Info.xPictureSize./results.Info.xPixel);
        sx = results.Info.xPixel;
        sy = results.Info.yPixel;
        
        [mx, my, wx, wy] = Gauss2D([],[],squeeze(sum(sum(imraw,3),4)),0);
        maxx = size(imraw);
        a(1)=round(mx-wx); a(2)=round(mx+wx); b(1)=round(my-wy); b(2)=round(my+wy);
        if a(1)<=0
            a(1)=1;
        end
        if a(2)>maxx(1)
            a(2)=maxx(1);
        end
        if b(1)<=0
            b(1)=1;
        end
        if b(2)>maxx(2)
            b(2)=maxx(2);
        end;
        
        imraw = imraw(a(1):a(2),b(1):b(2),:,:);
        [sx,sy,nx,ny] = size(imraw);
        mx = (nx-1)/2;
        my = (ny-1)/2;
        [x,y] = meshgrid(-my:my,-mx:mx);
        rx = (sx-1)/2;
        ry = (sy-1)/2;
        xx = -rx:rx;
        yy = -ry:ry;
        
        shift_size = 0.5*pixel_size/step_size;
        im = zeros(nx,ny);
        for jx=1:sx
            for jy = 1:sy
                im = im + interp2(x,y,squeeze(imraw(jx,jy,:,:)),x+xx(jx)*shift_size,y-yy(jy)*shift_size,'cubic',0);
            end
        end
        
        bx = ceil(shift_size*rx)+1;
        by = ceil(shift_size*ry)+1;
        im0 = squeeze(sum(sum(imraw(:,:,bx:end-bx+1,by:end-by+1),1),2));
        im0 = im0-min(im0(:));
        im = im(bx:end-bx+1,by:end-by+1);
        im = im-min(im(:));
        mim(cat(3,im0,im))
        eval(['print -dpng -r300 ''' dname fname(j).name(1:end-4) '_ScanImage'''])
    end
else
    load([dname fname])
    imraw = reshape(results.SHR,round([results.Info.xCCD results.Info.yCCD results.Info.xPixel results.Info.yPixel]));
    
    mim(sum(sum(imraw,3),4))
    eval(['print -dpng -r300 ''' dname fname(1:end-4) '_CCD'''])

    step_size = (results.Info.xPictureSize./results.Info.xPixel);
    sx = results.Info.xPixel;
    sy = results.Info.yPixel;
    
    [mx, my, wx, wy] = Gauss2D([],[],squeeze(sum(sum(imraw,3),4)),0);
    maxx = size(imraw);
    a(1)=round(mx-wx); a(2)=round(mx+wx); b(1)=round(my-wy); b(2)=round(my+wy);
    if a(1)<=0
        a(1)=1;
    end
    if a(2)>maxx(1)
        a(2)=maxx(1);
    end
    if b(1)<=0
        b(1)=1;
    end
    if b(2)>maxx(2)
        b(2)=maxx(2);
    end;
    
    imraw = imraw(a(1):a(2),b(1):b(2),:,:);
    [sx,sy,nx,ny] = size(imraw);
    mx = (nx-1)/2;
    my = (ny-1)/2;
    [x,y] = meshgrid(-my:my,-mx:mx);
    rx = (sx-1)/2;
    ry = (sy-1)/2;
    xx = -rx:rx;
    yy = -ry:ry;

    shift_size = 0.5*pixel_size/step_size;
    im = zeros(nx,ny);
    for jx=1:sx
        for jy = 1:sy
            im = im + interp2(x,y,squeeze(imraw(jx,jy,:,:)),x+xx(jx)*shift_size,y-yy(jy)*shift_size,'cubic',0);
        end
    end
    
    bx = ceil(shift_size*rx)+1;
    by = ceil(shift_size*ry)+1;
    im0 = squeeze(sum(sum(imraw(:,:,bx:end-bx+1,by:end-by+1),1),2));
    im0 = im0-min(im0(:));
    im = im(bx:end-bx+1,by:end-by+1);
    im = im-min(im(:));
    mim(cat(3,im0,im))
    eval(['print -dpng -r300 ''' dname fname(1:end-4) '_ScanImage'''])
end

if 0

    [nx, ny] = size(im);
    mx = (nx-1)/2;
    my = (ny-1)/2;
    [x,y] = meshgrid(-my:my,-mx:mx);
    
    % [a,b] = PickPixel;
    % clf;
    % plot(1:ny-2*by+1,im(bx+a,by:end-by),1:ny-2*by+1,res(bx+a,by:end-by,ind))
    
    clear res qfac
    NAv = 1.2;
    for k=1:length(NAv)
        NA = NAv(k);
        n0 = 1.33;
        n = 1.33;
        n1 = 1.33;
        d0 = [];
        d = 0;
        d1 = [];
        lamex = 0.63;
        lamem = 0.67;
        fd = 3e3;
        over = 7e3;
        pixel = 24;
        mag = pixel/step_size;
        pic = 1;
        focus = 0;
        zpos = 0;
        if mod(nx,2)==0
            mx = min([mx 80.5]);
        else
            mx = min([mx 80]);
        end
        if mod(ny,2)==0
            my = min([my 80.5]);
        else
            my = min([my 80]);
        end
        [mdf,exc] = MDFWideFieldMicroscope(NA,n0,n,n1,d0,d,d1,lamem,mag/2,focus,pixel,zpos,[mx,my],lamex,fd,over);
        [mdf,exc0] = MDFWideFieldMicroscope(NA,n0,n,n1,d0,d,d1,lamem,mag,focus,pixel,zpos,[mx,my],lamex,fd,over);
        exc = [zeros((nx-2*mx-1)/2,size(exc,2)); exc; zeros((nx-2*mx-1)/2,size(exc,2))];
        exc = [zeros(size(exc,1),(ny-2*my-1)/2), exc, zeros(size(exc,1),(ny-2*my-1)/2)];
        exc0 = [zeros((nx-2*mx-1)/2,size(exc0,2)); exc0; zeros((nx-2*mx-1)/2,size(exc0,2))];
        exc0 = [zeros(size(exc0,1),(ny-2*my-1)/2), exc0, zeros(size(exc0,1),(ny-2*my-1)/2)];
        mdf = [zeros((nx-2*mx-1)/2,size(mdf,2)); mdf; zeros((nx-2*mx-1)/2,size(mdf,2))];
        mdf = [zeros(size(mdf,1),(ny-2*my-1)/2), mdf, zeros(size(mdf,1),(ny-2*my-1)/2)];
        mx = (nx-1)/2;
        my = (ny-1)/2;
        otm1 = abs(fftshift(fft2(exc)));
        otm1 = otm1/max(otm1(:));
        %     otm2 = abs(fftshift(fft2(mdf)));
        %     otm2 = interp2((-2*my:2:2*my)',-2*mx:2:2*mx,otm2,(-my:my)',-mx:mx);
        %     otm2 = otm2/max(otm2(:));
        otm3 = abs(fftshift(fft2(mConv2(exc0,mdf))));
        otm3 = interp2((-2*my:2:2*my)',-2*mx:2:2*mx,otm3,(-my:my)',-mx:mx);
        otm3 = otm3/max(otm3(:));
        
        %weight = 1./(0.01+otm2).*(otm2>0.01);
        weight = otm1./(0.01+otm3).*(otm1>0.01);
        res(:,:,k) = abs(ifft2(ifftshift(fftshift(fft2(im)).*weight)));
        
        tmp = abs(fftshift(fft2(res(:,:,k))));
        qfac(k) = sum(sum((x.^2+y.^2).*tmp))/sum(sum(tmp));
        
        mim(cat(3,im,res(:,:,k)))
    end
    s = 1:length(NAv);
    ind = s(qfac==min(qfac));
    
    clf
    bx = 10; by = 10;
    if nx>ny
        [mx,my] = size(im0(bx:end-bx+1,by:end-by+1));
        tmp = [im0(bx:end-bx+1,by:end-by+1)/max(max(im0)) im(bx:end-bx+1,by:end-by+1)/max(max(im)) squeeze(res(bx:end-bx+1,by:end-by+1,ind))/max(max(res(:,:,ind)))];
        mim(tmp)
        colormap jet
        hold on
        plot([1 1]*my+1,[0 1]*(mx+0.5),'y')
        plot([2 2]*my+1,[0 1]*(mx+0.5),'y')
        plot([5 5+1/step_size],[5 5],'y','linewidth',3)
        hold off
        text((size(im,2)-2*by+2)/2,-size(im,1)/10,'raw image','HorizontalAlignment','center')
        text(3*(size(im,2)-2*by+2)/2,-size(im,1)/10,'ISM image','HorizontalAlignment','center')
        text(5*(size(im,2)-2*by+2)/2,-size(im,1)/10,{'Fourier filtered','ISM image'},'HorizontalAlignment','center')
    else
        [mx,my] = size(im0(bx:end-bx+1,by:end-by+1)');
        tmp = [im0(bx:end-bx+1,by:end-by+1)/max(max(im0)) im(bx:end-bx+1,by:end-by+1)/max(max(im)) squeeze(res(bx:end-bx+1,by:end-by+1,ind))/max(max(res(:,:,ind)))];
        mim(tmp)
        colormap jet
        hold on
        plot([1 1]*mx+0.5,[0 1]*(my+1),'y')
        plot([2 2]*mx+0.5,[0 1]*(my+1),'y')
        plot([5 5+1/step_size],[5 5],'y','linewidth',3)
        hold off
        text((size(im,1)-2*bx+2)/2,-size(im,2)/10,'raw image','HorizontalAlignment','center')
        text(3*(size(im,1)-2*bx+2)/2,-size(im,2)/10,'ISM image','HorizontalAlignment','center')
        text(5*(size(im,1)-2*bx+2)/2,-size(im,2)/10,{'apodized ISM image'},'HorizontalAlignment','center')
    end
    
end


