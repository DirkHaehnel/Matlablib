% Frap-Analysis

im1 = double(imread('E:\Daten\Iris\050413\1FRAP.tif',2));
mim(flipud(im1))
[a,b] = rubberband(gca,'-mode','line');
pos = ones(101,1)*a + (0:100)'/100*(b-a);
p1 = sqrt(sum((pos-ones(101,1)*pos(1,:)).^2,2));
z1 = interp2(1:size(im1,2),1:size(im1,1),im1,pos(:,1),pos(:,2),'cubic');
im2 = double(imread('E:\Daten\Iris\050413\1FRAP.tif',47));
mim(flipud(im2))
z2 = interp2(1:size(im2,2),1:size(im2,1),im2,pos(:,1),pos(:,2),'cubic');

pos = ginput(4);
ratio = (pos(4,2)-pos(1,2))/(pos(3,2)-pos(1,2));
rad = abs(pos(2,1)-pos(4,1));
tt = interp1(y(1,:),t,ratio,'cubic');


return

clear z y
x = 0:0.01:2;
t = 0.02:0.02:1;
for j=1:length(t)
    [ytmp,ztmp] = Frap(x,t(j));
    y(:,j) = ytmp';
    z(:,j) = ztmp';
end
save FrapTheory x t y z