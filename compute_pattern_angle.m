z = 0;
NA = 1.49;
n0 = 1.52;
n = 1.52;
n1 = 1.46;
d0 = [];
d = 0.01;
d1 = [];
lamem = 0.69;
mag = 100;
focus = 0:0.1:2;
atf = [];
ring = [];
pixel = 8;
pic = 1;
nn = [50 50];
al=pi/2;
be=0;

for i=1:length(focus)
    [im, mask] = PatternGeneration_Angle(z, NA, n0, n, n1, d0, d, d1, lamem, mag, focus(i), atf, ring, pixel, nn, be, al, pic);
      image{i}=im;
      maskim{i}=mask;
      clear mask
      clear im
end

% c=0
% for i=1:ceil(sqrt(focus))
%     for j=1:round(sqrt(focus))
%         c=c+1;
%         subplot(ceil(sqrt(length(focus))), round(sqrt(length(focus))),c)
%         
%     end
% end


for i=1:length(focus)
    maxi(i)=max(max(image{i}));
    image_normal{i}=image{i}./maxi(i);
end
 image_combine=cell2mat(image_normal);
 
 
 mim(image_combine);
