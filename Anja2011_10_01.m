dirname = 'm:\SOFI\FASL_Tests\Images\20110922_2\SOFI\';
dirname = 'm:\SOFI\FASL_Tests\Images\20110923_2\SOFI\';
fnames = dir([dirname '*.tif']);
nn = length(fnames);

for j=1:nn/2
    imr(:,:,j) = double(imread([dirname fnames(2*j-1).name]));
    imb(:,:,j) = double(imread([dirname fnames(2*j).name]));
end

stepxy = 0.100/2;
stepz = 0.150;
[x y z] = meshgrid(stepxy*(1:size(imr,1)),stepxy*(1:size(imr,2)),stepz*(1:size(imr,3)));

ind = 0:255;
for j=1:nn/2
    tst = imr(:,:,j);
    v = mHist(tst(:),ind);
    tst = imr(:,:,j) - ind(v==max(v));
    tst(tst<0) = 0;
    tst = tst/max(tst(:))*255;
    imr2(:,:,j) = tst;
    tst = imb(:,:,j);
    v = mHist(tst(:),ind);
    tst = imb(:,:,j) - ind(v==max(v));
    tst(tst<0) = 0;
    tst = tst/max(tst(:))*255;
    imb2(:,:,j) = tst;
end

if 0
    thr = 70;
    thb = 30;
    close all
    pr=patch(isosurface(x,y,z,imr2,thr));
    isonormals(x,y,z,imr2,pr);
    set(pr,'FaceColor','red','EdgeColor','none','FaceAlpha',0.8);
    pb=patch(isosurface(x,y,z,imb2,thb));
    isonormals(x,y,z,imb2,pb);
    set(pb,'FaceColor','blue','EdgeColor','none','FaceAlpha',0.8);
    daspect([1 1 1])
    view(3); axis vis3d; axis off
    line(x(1,:,1),y(1,:,1),z(1,:,1),'color',[0.7 0.7 0.7],'linewidth',1)
    line(x(1,:,1),y(end,:,1),z(1,:,1),'color',[0.7 0.7 0.7],'linewidth',1)
    line(x(1,:,1),y(1,:,1),z(1,:,end),'color',[0.7 0.7 0.7],'linewidth',1)
    line(x(1,:,1),y(end,:,1),z(1,:,end),'color',[0.7 0.7 0.7],'linewidth',1)
    line(x(:,1,1),y(:,1,1),z(:,1,1),'color',[0.7 0.7 0.7],'linewidth',1)
    line(x(:,end,1),y(:,1,1),z(:,1,1),'color',[0.7 0.7 0.7],'linewidth',1)
    line(x(:,1,1),y(:,1,1),z(:,1,end),'color',[0.7 0.7 0.7],'linewidth',1)
    line(x(:,end,1),y(:,1,1),z(:,1,end),'color',[0.7 0.7 0.7],'linewidth',1)
    line(squeeze(x(1,1,:)),squeeze(y(1,1,:)),squeeze(z(1,1,:)),'color',[0.7 0.7 0.7],'linewidth',1)
    line(squeeze(x(1,end,:)),squeeze(y(1,1,:)),squeeze(z(1,1,:)),'color',[0.7 0.7 0.7],'linewidth',1)
    line(squeeze(x(1,1,:)),squeeze(y(end,1,:)),squeeze(z(1,1,:)),'color',[0.7 0.7 0.7],'linewidth',1)
    line(squeeze(x(1,end,:)),squeeze(y(end,1,:)),squeeze(z(1,1,:)),'color',[0.7 0.7 0.7],'linewidth',1)
    camlight
    lighting gouraud
    
    for j=1:36
        view([10*j,25])
        %view([cos(j*pi/18),-cos(j*pi/18),sin(j*pi/18)])
        eval(['print -dpng -r300 neuron1_' mint2str(j,2)])
    end
end

if 0
    for j=2:nn/2
        thr = 70;
        thb = 30;
        close all
        pr=patch(isosurface(x(:,:,1:j),y(:,:,1:j),z(:,:,1:j),imr2(:,:,1:j),thr));
        isonormals(x(:,:,1:j),y(:,:,1:j),z(:,:,1:j),imr2(:,:,1:j),pr);
        set(pr,'FaceColor','red','EdgeColor','none','FaceAlpha',0.8);
        pb=patch(isosurface(x(:,:,1:j),y(:,:,1:j),z(:,:,1:j),imb2(:,:,1:j),thb));
        isonormals(x(:,:,1:j),y(:,:,1:j),z(:,:,1:j),imb2(:,:,1:j),pb);
        set(pb,'FaceColor','blue','EdgeColor','none','FaceAlpha',0.8);
        patch(x(:,:,1),y(:,:,1),z(:,:,1),'FaceColor','w','EdgeColor','none')
        axis image
        axis tight
        whitebg('w')
        box on
        camlight
        lighting gouraud
        eval(['print -dpng -r300 neuron1_' mint2str(j-1,2)])
        eval(['print -dpng -r300 neuron1_' mint2str(nn-j,2)])
    end
end

