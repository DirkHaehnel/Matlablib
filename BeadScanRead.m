% BeadScanRead

clear mx1 my1 wx1 wy1 amp1 mx2 my2 wx2 wy2 amp2 d zd im1 im2

path = 'D:\Joerg\Doc\Fcs\2Focus\Anastasia\scan\';

names = dir([path 'Atto fcs sur d20 11.t3r']);

read_fun = 'SCX200Read';
%read_fun = 'PIE710Read';

close all hidden
pix = 200e3/2^12;
for j=1:length(names)
    %zd(j) = str2num(names(j).name(5:7));
    zd = 0;
end
[zd ind] = sort(zd);
names = names(ind);
for j=1:length(names)
    fname = [path names(j).name];

    eval(['[tag, tim, tch, bin] = ' read_fun '(fname);'])

    if j==1
        semilogy(bin,tch);
        pos = ginput(2);
    end

    filter = bin>pos(1,1) & bin<pos(2,1);
    eval(['[tag1, tim1] = ' read_fun '(fname,[filter'' filter'']);'])
    filter = bin<pos(1,1) | bin>pos(2,1);
    eval(['[tag2, tim2] = ' read_fun '(fname,[filter'' filter'']);'])

    if j==1
        im1(:,:,j) = sum(tim1,3);
        im2(:,:,j) = sum(tim2,3);
    else
        mx = size(tim1,1);
        my = size(tim1,2);
        im1(1:mx,1:my,j) = sum(tim1,3);
        im2(1:mx,1:my,j) = sum(tim2,3);
    end
end

return

for j=1:1
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

return

al=1.04; be=-200; 
plot(al*exa.z,exa.wx1,al*exa.z,exa.wy1,exb.z/al-be,exb.wx1,exb.z/al-be,exb.wy1)
plot(al*exa.z,exa.wx2,al*exa.z,exa.wy2,exb.z/al-be,exb.wx2,exb.z/al-be,exb.wy2)
plot(al*exa.z,exa.amp1,exb.z/al-be,exb.amp1)
plot(al*exa.z,exa.amp2,exb.z/al-be,exb.amp2)
plot(al*exa.z,exa.amp1,exb.z/al-be,exb.amp1,al*exa.z,exa.amp2,exb.z/al-be,exb.amp2)
plot(al*exa.z,sqrt((exa.mx1-exa.mx2).^2+(exa.my1-exa.my2).^2),exb.z/al-be,sqrt((exb.mx1-exb.mx2).^2+(exb.my1-exb.my2).^2))

return

w0=0.4;z0=2.7; plot(zv(ind),wx1(ind),zv(ind),wy1(ind),zv(ind),w0*sqrt(1+(zv(ind)-z0).^2*0.65/pi/w0^2)*1e3);
plot(zv,(ampx1+ampy1).*wx1.*wy1/max((ampx1+ampy1).*wx1.*wy1),zv,(ampx2+ampy2).*wx2.*wy2/max((ampx2+ampy2).*wx2.*wy2))

names = dir('D:\Daten\Dertinger\102605PSFScan\z_e-5_p61_*.t3r');
for j=1:length(names)
    fname = ['D:\Daten\Dertinger\102605PSFScan\' names(j).name];
    [bla, bla, bla, head(j)] = tttrread(fname,[1 0]);
end

joerg = {'psf52', 'psf56', 'psf61', 'psf65'};
amp = [];
for jj=1:length(joerg)
    load(joerg{jj})
    amp = [amp ampx1+ampy1 ampx2+ampy2];
end