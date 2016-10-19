namefl = 'c:\Joerg\Doc\Ales\Data2007-10-13\B10_f.txt';
namein = 'c:\Joerg\Doc\Ales\Data2007-10-13\B10_e.txt';

%t = [0,0,1,1,1,1,2,2,2,5,5,5,8,8,10];
t = [0,1,1,1,2,2,4,4,7,10];
nx = 60;
ny = 60;

fin = fopen(namefl,'r');
y = fscanf(fin,'%f ',inf);
fclose(fin);
fin = fopen(namein,'r');
x = fscanf(fin,'%f ',inf);
fclose(fin);

x = reshape(x,nx,length(t),ny);
y = reshape(y,nx,length(t),ny);

[b,a]=meshgrid(1:size(y,1),1:size(y,3));
a = a(y(:,end,:)==max(y(:)));
b = b(y(:,end,:)==max(y(:)));

%plot(x(a,:,b),y(a,:,b))

plot(mean(mean(x,3),1),mean(mean(y,3),1))
xlabel('rel. excitation intensity')
ylabel('rel. fluorescence');

return

plot(mean(mean(x,3),1)/max(mean(mean(x,3),1)),[mean(mean(y(45:55,:,45:55),3),1)/max(mean(mean(y(45:55,:,45:55),3),1));mean(mean(y,3),1)/max(mean(mean(y,3),1))])

close; [mx, my, wx0, wy0] = Gauss2D(squeeze(y(a+(-15:15),end,b+(-15:15))./x(a+(-15:15),end,b+(-15:15))));
close; [mx, my, wx1, wy1] = Gauss2D(squeeze(mean(y(a+(-15:15),t==4,b+(-15:15))./x(a+(-15:15),t==4,b+(-15:15)),2))-squeeze(y(a+(-15:15),end,b+(-15:15))./x(a+(-15:15),end,b+(-15:15))));