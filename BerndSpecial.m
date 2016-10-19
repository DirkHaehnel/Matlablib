% program for OTF determination

close all

%load('m:\ISM\2010-02\2010-02-05\TS100_Stack-10.mat')

a = 74:103;
b = 130:158;
% a = 136:156;
% b = 114:133;

aa = 12:20;
bb = 8:16;

imraw = results.SHR;
tst = squeeze(sum(sum(imraw,3),4)); mim(tst)
[ind,ind] = max(tst(:));
clust = mCluster(tst>mean(tst(:))); 
clust = clust==clust(ind);
[sx,sy,nx,ny] = size(imraw);
mx = (nx-1)/2;
my = (ny-1)/2;
[x,y] = meshgrid(-my:my,-mx:mx);
rx = (sx-1)/2;
ry = (sy-1)/2;
[xx,yy] = meshgrid(-rx:rx,-ry:ry);
ym = xx(ind)+rx+1;
xm = yy(ind)+ry+1;
xx = xx - yy(ind);
yy = yy - xx(ind);

pixel_size = 0.08;
step_size = 1e-3*results.Info.xPictureSize/results.Info.xPixel;

close all; 
for j=1:length(aa) 
    for k=1:length(bb) 
        [mx(j,k) my(j,k) wx(j,k) wy(j,k)] = Gauss2D(squeeze(imraw(aa(j),bb(k),a,b))); 
        drawnow; 
    end; 
end
sx = mx - mean(mx(:));
sy = my - mean(my(:));

xx = xx(1,:);
yy = yy(:,1)';
im = zeros(nx,ny);
im0 = im;
mm = mean(mean(squeeze(sum(sum(imraw,3),4))));
for jx=aa(1):aa(end)
    for jy = bb(1):bb(end)
        im = im + interp2(x,y,squeeze(imraw(jx,jy,:,:)),x+sx(jx-aa(1)+1,jy-bb(1)+1),y+sy(jx-aa(1)+1,jy-bb(1)+1),'cubic',0);
        im0 = im0 + squeeze(imraw(jx,jy,:,:));
    end
end

bx = ceil(abs(shift_size_x)*rx)+1;
by = ceil(abs(shift_size_y)*ry)+1;
im0 = im0(bx:end-bx+1,by:end-by+1);
im0 = im0-min(im0(:));
im = im(bx:end-bx+1,by:end-by+1);
im = im-min(im(:));
mim(cat(3,im0,im))
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
    lamex = 0.64;
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
im0 = im0(bx:end-bx+1,by:end-by+1);
im = im(bx:end-bx+1,by:end-by+1);
res = squeeze(res(bx:end-bx+1,by:end-by+1,ind));
%if nx>ny
    [mx,my] = size(im0);
    tmp = [im0/max(max(im0)) im/max(max(im)) res/max(max(res))];
    mim(tmp)
    %colormap jet
    hold on
    plot([1 1]*mx+0.5,[0 1]*(my+1),'y')
    plot([2 2]*mx+0.5,[0 1]*(my+1),'y')
    plot([5 5+1/step_size],[5 5],'y','linewidth',3)
    hold off
    text((size(im,2)+2)/2,-size(im,1)/10,'raw image','HorizontalAlignment','center')
    text(3*(size(im,2)+2)/2,-size(im,1)/10,'ISM image','HorizontalAlignment','center')
    text(5*(size(im,2)+2)/2,-size(im,1)/10,{'Fourier filtered','ISM image'},'HorizontalAlignment','center')
% else
%     [mx,my] = size(im0');
%     tmp = [im0'/max(max(im0)) im'/max(max(im)) res'/max(max(res))];
%     mim(tmp)
%     %colormap jet
%     hold on
%     plot([1 1]*mx+0.5,[0 1]*(my+1),'y')
%     plot([2 2]*mx+0.5,[0 1]*(my+1),'y')
%     plot([5 5+1/step_size],[5 5],'y','linewidth',3)
%     hold off
%     text((size(im,1)+2)/2,-size(im,2)/10,'raw image','HorizontalAlignment','center')
%     text(3*(size(im,1)+2)/2,-size(im,2)/10,'ISM image','HorizontalAlignment','center')
%     text(5*(size(im,1)+2)/2,-size(im,2)/10,{'apodized ISM image'},'HorizontalAlignment','center')
% end

return

im=squeeze(sum(sum(imraw,1),2)); for j=1:size(imraw,3) for k=1:size(imraw,4) subplot(121); mim(im); hold on; plot(k,j,'oc'); hold off; subplot(122); mim(imraw(:,:,j,k)'); pause(0.01); end; end

return

[a,b]=PickPixel
tmp = [im0(bx+b+(-20:20),by+a+(-20:20))'/max(max(im0)) im(bx+b+(-20:20),by+a+(-20:20))'/max(max(im)) squeeze(res(bx+b+(-20:20),by+a+(-20:20),ind))'/max(max(res(:,:,ind)))];
mim(tmp)
colormap jet
hold on
plot([1 1]*size(tmp,2)/3,[0 1]*size(tmp,1),'y')
plot([2 2]*size(tmp,2)/3,[0 1]*size(tmp,1),'y')
plot([5 5+1/step_size],[2 2],'y','linewidth',3)
hold off
text(size(tmp,2)/3/2,-size(im,2)/15,'raw image','HorizontalAlignment','center')
text(3*size(tmp,2)/3/2,-size(im,2)/15,'ISM image','HorizontalAlignment','center')
text(5*size(tmp,2)/3/2,-size(im,2)/15,{'apodized ISM image'},'HorizontalAlignment','center')

return

mim(im)
[b,a] = ginput(2);
a=round(a);b=round(b);
[mx1,my1,wx1,wy1]=GaussEllipse([],[],im(a(1):a(2),b(1):b(2)));
for k=1:length(NAv) [mx2(k),my2(k),wx2(k),wy2(k)]=Gauss2D(res(a(1):a(2),b(1):b(2),k)); end
plot(NAv,wx2,NAv,wy2,NAv,ones(size(NAv))*wx1,'ro',NAv,ones(size(NAv))*wy1,'bo',NAv,ones(size(NAv))*wx1/sqrt(2),'r--',NAv,ones(size(NAv))*wy1/sqrt(2),'b--')
ax = axis;
hold on
plot([1 1]*NAv(ind),ax(3:4),':k');
hold off

return

[sx,sy,nx,ny] = size(imraw);
tst=zeros(sx*nx,sy*ny);
for kx=1:sx for ky=1:sy tst((kx-1)*nx+1:kx*nx,(ky-1)*ny+1:ky*ny)=imraw(kx,ky,:,:); end; end

tst=zeros(sx*nx,sy*ny);
for kx=1:sx for ky=1:sy tst((kx-1)*nx+1:kx*nx,(ky-1)*ny+1:ky*ny)=im(kx,ky,:,:)/mean(mean(im(kx,ky,:,:))); end; end

tst=zeros(sx*nx,sy*ny);
for kx=1:sx for ky=1:sy tst((kx-1)*nx+1:kx*nx,(ky-1)*ny+1:ky*ny)=abs(fftshift(fft2(squeeze(im(kx,ky,:,:))-mean(squeeze(im(kx,ky,:)))))); end; end

return

im = zeros(nx,ny,10);
for k=1:10
    shift_size = k/10*pixel_size/step_size;
    for jx=1:sx
        for jy = 1:sy
            im(:,:,k) = im(:,:,k) + interp2(x,y,squeeze(imraw(jx,jy,:,:)),x-(xx(1,jx)+3)*shift_size,y+yy(jy,1)*shift_size,'cubic',0);
        end
    end
end

return

% z_scan
for j=1:length(fnames) load(['m:\ISM\2010-02\2010-02-09\Neuer Ordner\' fnames(j).name]); BerndFourier; final(j).im0=im0; final(j).im=im; final(j).res=res; final(j).shift_size_x=shift_size_x; final(j).shift_size_y=shift_size_y; end
close all
while 1==1
    for j=1:16
        tmp = [final(j).im0/max(max(final(j).im0)) final(j).im/max(max(final(j).im)) final(j).res/max(max(final(j).res))];
        mim(tmp)
        hold on
        plot([1 1]*mx+0.5,[0 1]*(my+1),'y')
        plot([2 2]*mx+0.5,[0 1]*(my+1),'y')
        plot([5 5+1/step_size],[5 5],'y','linewidth',3)
        hold off
        drawnow;
    end
end
