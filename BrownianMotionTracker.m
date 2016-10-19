function [pthx, pthy, mask] = BrownianMotionTracker(z,tsh,mask)

if nargin<2 || isempty(tsh) 
    tsh = 2;
end

if nargin<3 || isempty(mask)
    mim(z(:,:,1));
    [a,b] = ginput(2);
    a = round(a); b = round(b);
    mask = z(b(1):b(2),a(1):a(2));
    [mx,my,wx,wy] = Gauss2D(mask);
    mm = ceil(max(size(mask))/2);
    [mx,my] = meshgrid(-mm:mm,-mm:mm);
    mask = exp(-2*mx.^2/wx.^2-2*my.^2/wy.^2);
end
    
[err,bim,cim,sim,xc,yc] = FindPattern(z(:,:,1),mask,[],[],[],tsh);
a(1:length(xc),1) = xc';
b(1:length(xc),1) = yc';
mim(z(:,:,1)); hold on; plot(xc,yc,'o'); hold off; drawnow
for k=2:size(z,3)
    [err,bim,cim,sim,xc,yc] = FindPattern(z(:,:,k),mask,[],[],[],tsh);
    a(1:length(xc),k) = xc';
    b(1:length(xc),k) = yc';
    mim(z(:,:,k)); hold on; plot(xc,yc,'o'); hold off; drawnow
end
    
[pthx,pthy] = Neighbor(a./(a>0),b./(b>0));
maxdist = 10;
len = length(pthx);
cnt = len;
for j=1:len 
    ind = 1:length(pthx{j});
    ind = ind(abs(diff(pthx{j}))>maxdist | abs(diff(pthy{j}))>maxdist);
    if ~isempty(ind)
        ind = [ind length(pthx{j})];
        for k=1:length(ind)-1
            cnt = cnt+1;    
            pthx{cnt} = pthx{j}(ind(k)+1:ind(k+1));
            pthy{cnt} = pthy{j}(ind(k)+1:ind(k+1));
        end
        pthx{j} = pthx{j}(1:ind(1));
        pthy{j} = pthy{j}(1:ind(1));
    end
end
for j=1:length(pthx)
    len(j) = length(pthx{j});
end
pthx(len==1)=[];
pthy(len==1)=[];
len(len==1)=[];
[ind,ind]=sort(len);
pthx=pthx(ind);
pthy=pthy(ind);

frb = ['mcrgb'];
plot(pthx{1},pthy{1},frb(1)); 
axis([0 size(z,2) 0 size(z,1)]); 
hold on; 
for j=2:length(pthx) 
    plot(pthx{j},pthy{j},frb(mod(j-1,5)+1)); 
end; 
hold off
axis image

%  for j=1:length(pthx) pthx{j}=pthx{j}-pthx{j}(1); pthy{j}=pthy{j}-pthy{j}(1); end