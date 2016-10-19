% CoverSlideAnalysis

clear all

path = 'D:\Joerg\Doc\Fcs\2Focus\Coverslide\2006-02-21\';

names = dir([path 'C_R_*b.mat']);

close all hidden
for j=1:length(names)
    zd(j) = str2num(names(j).name(5:end-5));
end
[zd ind] = sort(zd);
names = names(ind);
for j=1:length(names)
    fname = [path names(j).name];
    load(fname);
    h(j) = head;
    para(j) = parameters;
    auto(:,:,j) = [res.auto(:,1,2)+res.auto(:,2,1),res.auto(:,3,4)+res.auto(:,4,3),res.auto(:,1,4)+res.auto(:,4,1)+res.auto(:,2,3)+res.auto(:,3,2)];
end
t = res.autotime;

for j=1:size(auto,3)
    p1(:,j) = simplex('rigler',[1e-3 1e-2],[0 0],[],[],[],t,auto(:,1,j),[],1);
    for k=1:2
        p1(:,j) = simplex('rigler',p1(:,j),[0 0],[],[],[],t,auto(:,1,j),[],1);
    end
    p2(:,j) = simplex('rigler',[1e-3 1e-2],[0 0],[],[],[],t,auto(:,2,j),[],1);    
    for k=1:2
        p2(:,j) = simplex('rigler',p2(:,j),[0 0],[],[],[],t,auto(:,2,j),[],1);
    end
    px(:,j) = simplex('rigler',[1e-3 1e-2],[0 0],[],[],[],t,auto(:,3,j),[],1);        
    for k=1:2
        px(:,j) = simplex('rigler',px(:,j),[0 0],[],[],[],t,auto(:,3,j),[],1);
    end
end

plot(zd,p1(1,:)+p2(1,:),zd,polyval(polyfit(zd,p1(1,:)+p2(1,:),3),zd))
xlabel('cover slide adjustment collar [\mum]')
ylabel('diffusion time [s]');
 
tmp = polyval(polyfit(zd,1./(p1(1,:)+p2(1,:)),3),zd);
plot(zd,1./(p1(1,:)+p2(1,:))/max(tmp),zd,tmp/max(tmp))
xlabel('cover slide adjustment collar [\mum]')
ylabel('rel. diffusion coefficient');

plot(zd,min(p1(1,:))./p1(1,:),zd,min(p2(1,:))./p2(1,:))
plot(zd,0.41e-4^2./2./(px(1,:)-(p1(1,:)+p2(1,:))/2),'o-')

for j=1:size(auto,3)
    pp(:,j) = simplex('Rigler2f',[1e-3 1e-2 1e-2],[0 0 0],[],[],[],t,squeeze(auto(:,1:2,j)),auto(:,3,j)/2);
    for k=1:2
        pp(:,j) = simplex('Rigler2f',pp(:,j),[0 0 0],[],[],[],t,squeeze(auto(:,1:2,j)),auto(:,3,j)/2);
    end
end
plot(zd,0.41e-4^2./pp(3,:)/4,'o-')
xlabel('cover slide adjustment collar [\mum]')
ylabel('diffusion coefficient [cm^2/s]');
