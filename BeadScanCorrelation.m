clear mx1 my1 wx1 wy1 amp1 mx2 my2 wx2 wy2 amp2 d zd im1 im2

path = 'D:\Joerg\Doc\Fcs\2Focus\PSF\PSF_Iris_2006-08-17\';

names = dir([path 'z*a.t3r']);

close all hidden
pix = 200e3/2^12;
for j=1:length(names)
    zd(j) = str2num(names(j).name(3:5));
    %zd = 0;
end
[zd ind] = sort(zd);
names = names(ind);
for j=1:length(names)
    fname = [path names(j).name];

    [tag, tim, tch, bin] = SCX200Read(fname);

    if j==1
        semilogy(bin,tch);
        pos = ginput(2);
    end

    filter = bin>pos(1,1) & bin<pos(2,1);
    [tag1, tim1] = SCX200Read(fname,[filter' filter']);
    filter = bin<pos(1,1) | bin>pos(2,1);
    [tag2, tim2] = SCX200Read(fname,[filter' filter']);

    if j==1
        im1(:,:,j) = sum(tim1,3);
        im2(:,:,j) = sum(tim2,3);
    else
        mx = size(tim1,1);
        my = size(tim1,2);
        im1(1:mx,1:my,j) = sum(tim1,3);
        im2(1:mx,1:my,j) = sum(tim2,3);
    end
    
    mim(cat(3,sum(tim1,3),sum(tim2,3)))
    subplot(121);
    [mx1(j),my1(j),wx1(j),wy1(j),amp1(j)] = Gauss2D([],[],sum(tim1,3),pi/4);
    subplot(122);
    [mx2(j),my2(j),wx2(j),wy2(j),amp2(j)] = Gauss2D([],[],sum(tim2,3),pi/4);

    if j>1
        if ~(my==size(im1,2))
            tmp = wx1(j); wx1(j) = wy1(j); wy1(j) = tmp;
            tmp = wx2(j); wx2(j) = wy2(j); wy2(j) = tmp;
        end
    end
    
    d(j) = sqrt((mx1(j)-mx2(j))^2+(my1(j)-my2(j))^2)*pix;

%         suptitle(['distance = ' num2str(d(j)) ' nm'])
%         eval(['print -dpng -zbuffer ' fname(1:end-4)])
end
ex.wx1 = pix*wx1;
ex.wx2 = pix*wx2;
ex.wy1 = pix*wy1;
ex.wy2 = pix*wy2;
ex.mx1 = pix*mx1;
ex.mx2 = pix*mx2;
ex.my1 = pix*my1;
ex.my2 = pix*my2;
ex.amp1 = amp1;
ex.amp2 = amp2;
ex.im1 = im1;
ex.im2 = im2;
ex.z = zd;
ex.d = d;

% correlation analysis

clear auto1 auto2 auto12 auto21
len = min([size(ex.im1,1) size(ex.im1,2)]);
for k=1:1:length(names) 
    for j=1:len-1 
        auto1(j,k,1) = sum(sum(ex.im1(1:len-j,:,k).*ex.im1(j+1:len,:,k),1)/(len-j),2); 
        auto2(j,k,1) = sum(sum(ex.im2(1:len-j,:,k).*ex.im2(j+1:len,:,k),1)/(len-j),2); 
        auto12(j,k,1) = sum(sum(ex.im1(1:len-j,:,k).*ex.im2(j+1:len,:,k),1)/(len-j),2); 
        auto21(j,k,1) = sum(sum(ex.im2(1:len-j,:,k).*ex.im1(j+1:len,:,k),1)/(len-j),2);         
        auto1(j,k,2) = sum(sum(ex.im1(:,1:len-j,k).*ex.im1(:,j+1:len,k),1)/(len-j),2); 
        auto2(j,k,2) = sum(sum(ex.im2(:,1:len-j,k).*ex.im2(:,j+1:len,k),1)/(len-j),2); 
        auto12(j,k,2) = sum(sum(ex.im1(:,1:len-j,k).*ex.im2(:,j+1:len,k),1)/(len-j),2); 
        auto21(j,k,2) = sum(sum(ex.im2(:,1:len-j,k).*ex.im1(:,j+1:len,k),1)/(len-j),2);         
    end; 
end

close all
clear p1 p2
t = 1:30;
for j = 1:length(names);
    tmp = Simplex('BeadScanCorrelationFun',[10 10],[0 -inf],[],[],[],t,[auto1(t,j,1),auto2(t,j,1),auto12(t,j,1),auto21(t,j,1)]);
    p1(j) = tmp(2);
    tmp = Simplex('BeadScanCorrelationFun',[10 10],[0 -inf],[],[],[],t,[auto1(t,j,2),auto2(t,j,2),auto12(t,j,2),auto21(t,j,2)]);
    p2(j) = tmp(2);
end

plot(1:length(names),sqrt(p1.^2+p2.^2)*pix,1:length(names),ex.d)