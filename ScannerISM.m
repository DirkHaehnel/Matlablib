% program for OTF determination

clear all 
close all




[filename, pathname] = uigetfile('*.mat', 'Pick an ISM-File','\\jesrv\WG-Data\ISM\*');
if isequal(filename,0)
    return
else
    load(fullfile(pathname, filename));
    filenameList = dir([pathname '*.mat']);
end

if isstruct(filenameList)
    for j=1:length(filenameList)

load([pathname filenameList(j).name])

im=results.SHR(1:1420800) 
imraw=reshape(im,16,16,74,75)     


mim(squeeze(sum(sum(imraw,3),4)));
eval(['print -dpng -r300 ''' pathname filenameList(j).name '_psf.png'''])
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
%[a,b] = PickPixel;

imraw = imraw(a(1):a(2),b(1):b(2),:,:);
[sx,sy,nx,ny] = size(imraw);
mx = (nx-1)/2;
my = (ny-1)/2;
[x,y] = meshgrid(-my:my,-mx:mx);
rx = (sx-1)/2;
ry = (sy-1)/2;
xx = -rx:rx;
yy = -ry:ry;

%pixel_size = 80;
%step_size = (results.Info.xPictureSize./results.Info.xPixel);

pixel_size = 0.08;
step_size = 0.08;
%step_size = 1e-3*results.Info.xPictureSize/results.Info.xPixel


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
[nx, ny] = size(im);
mx = (nx-1)/2;
my = (ny-1)/2;
[x,y] = meshgrid(-my:my,-mx:mx);

%  [a,b] = PickPixel;
%  clf;
%  plot(1:ny-2*by+1,im(bx+a,by:end-by),1:ny-2*by+1,res(bx+a,by:end-by,ind))

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
    pixel =24;
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
    eval(['print -dpng -r300 ''' pathname filenameList(j).name '_image.png'''])
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
    eval(['print -dpng -r300 ''' pathname filenameList(j).name '_image.png'''])
end
    end
end
%end %for loop multicolor

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

% Verschiebung des Fokus-Punktes
tast=size(results.SHR);
for i=1:tast(1)
    for j=1:tast(2)
        pix = mConv2(squeeze(results.SHR(i,j,:,:)),squeeze(results.SHR(round(tast(1)/2),round(tast(2)/2),:,:)));
        [tst(1,:) tst(2,:)]=max(pix);
        [test(1) test(2)]=max(tst(1,:));
        maximum(i,j,1)=test(2);
        maximum(i,j,2)=tst(2,test(2));
    end
end

return

%Focus-Finder
clear ergebniss;
clear ergebniss1;
clear ergebniss2;


for i=1:numofpix
    name = ['\\jesrv\WG-Data\ISM\2011-09\testdaten\test\Data_ISM_X-Y_Cam_0_file1' num2str(i) '.mat'];
    load(name);
    ergebniss1(i,:,:) = squeeze(sum(sum(results.SHR,1),2));
    ergebniss2(i,:,:) = squeeze(sum(sum(results.SHR,3),4));
    [mx my wx wy amp]=Gauss2D([],[],squeeze(ergebniss2(i,5:40,10:45)));
    abm=size(results.SHR);
    ergebniss(i,:) = [results.Info.zPos (mx+5)./abm(1) (my+10)./abm(2) wx*0.08 wy*0.08 amp];
end

subplot(2,3,1);
plot(ergebniss(:,1),[ergebniss(:,2) ergebniss(:,3)]);
legend({'mx','my'});
xlabel('z-Pos / [nm]');
ylabel('Position / [nm]');
subplot(2,3,2);
plot(ergebniss(:,1),[ergebniss(:,4) ergebniss(:,5)]);
legend({'wx','wy'});
xlabel('z-Pos / [nm]');
ylabel('1/e-Gauss, / [\mum]');
subplot(2,3,3);
plot(ergebniss(:,1),ergebniss(:,6));
legend({'amp'});
xlabel('z-Pos / [nm]');
ylabel('Counts');
drawnow
while 1==1
for i=1:numofpix
    subplot(2,3,4)
    mim(squeeze(ergebniss2(i,:,:)))
    subplot(2,3,5)
    mim(squeeze(ergebniss1(i,:,:)))
    drawnow
end
end