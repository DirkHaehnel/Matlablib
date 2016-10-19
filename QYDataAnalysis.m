% QYDataAnalysis

load D:\Joerg\Doc\Fcs\Cavity\scan.mat

tag = tag(:,:,1);
tim = tim(:,:,1);
[x,y] = meshgrid(1:size(tag,2),1:size(tag,1));
close; 
p=Simplex('RingFit',[1 80 50 100],[],[],[],[],tag)
for j=1:4
    p = Simplex('RingFit',p,[],[],[],[],tag)
end
r = sqrt((x-p(1)).^2+(y-p(2)).^2);
bin = 0.5:max(r(:));
int = mHist(r(:),bin,tag(:))./mHist(r(:),bin);
tau = mHist(r(:),bin,tim(:))./mHist(r(:),bin);
bin([1 end]) = [];
int([1 end]) = [];
tau([1 end]) = [];
plot(bin,tau,bin,(int/max(int)-min(int/max(int)))*max(tau))