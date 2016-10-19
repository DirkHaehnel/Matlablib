%p1 = load('f:\MSC-orientation-analysis\1kpa-radial-intensity-profiles.txt');
%p1 = load('f:\MSC-orientation-analysis\11kpa-radial-intensity-profiles.txt');
p1 = load('f:\MSC-orientation-analysis\34kpa-radial-intensity-profiles.txt');
ind = p1(2:end,1)<p1(1:end-1,1);
t = 1:size(p1,1);
t = [0 t(ind) size(p1,1)];
clear f1; for j=1:length(t)-1 f1{j}=p1(t(j)+1:t(j+1),:); end
for j=1:length(f1) 
    res1(:,j)=interp1(f1{j}(:,1)/max(f1{j}(:,1)),f1{j}(:,2),x,'cubic'); 
end
figure
plot(x,res1./(ones(size(res1,1),1)*res1(1,:)),x,mean(res1./(ones(size(res1,1),1)*res1(1,:)),2),'o')

