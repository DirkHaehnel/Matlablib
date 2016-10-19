function [pthx, pthy] = Pavlik(name,tsh)

if nargin<2 || isempty(tsh) 
    tsh = 10;
end

if nargin<1 || isempty(name)
    [name, path] = uigetfile('*.tif');
    name = name(1:end-4);
    ind = findstr(name,'_');
    files = dir([path name(1:ind-1) '*.tif']);
end

for k=1:length(files)
    num(k) = str2num(files(k).name(ind+1:end-4));
end
[num, ind] = sort(num);

k=1;
z = double(imread([path files(ind(k)).name]));
mim(z);
[a,b] = ginput(2);
a = round(a); b = round(b);
mask = z(b(1):b(2),a(1):a(2));
[err,cim,sim,xc,yc] = FindPattern(z,mask,[],0,10);
a(1:length(xc),k) = xc';
b(1:length(xc),k) = yc';
mim(z); hold on; plot(xc,yc,'o'); hold off; drawnow
for k=2:length(files)
    z = double(imread([path files(ind(k)).name]));
    [err,cim,sim,xc,yc] = FindPattern(z,mask,[],0,10);
    
    a(1:length(xc),k) = xc';
    b(1:length(xc),k) = yc';
    mim(z); hold on; plot(xc,yc,'o'); hold off; drawnow
end
    
ngbr = neighbor(a,b);
tmp = mhist(ngbr,0:size(ngbr,1));
tmp(1,:) = [];
for j1=1:size(ngbr,2)
    for j2=1:size(ngbr,1)
        if tmp(j2,j1)>1
            ngbr(ngbr(:,j1)==j2,j1) = 0;
        end
    end
end
pthx = {[]}; pthy = {[]}; cnt = 1;
for j1=1:size(ngbr,2)
    for j2=1:size(ngbr,1)
        if ngbr(j2,j1)>0
            nxt = ngbr(j2,j1);
            pthx{cnt}(1:2) = [a(j2,j1) a(nxt,j1+1)];
            pthy{cnt}(1:2) = [b(j2,j1) b(nxt,j1+1)];                
            ngbr(j2,j1) = 0;
            tst = 1; jj = j1+1;
            while tst>0 & jj<=size(ngbr,2)
                if ngbr(nxt,jj)>0
                    nxtnxt = ngbr(nxt,jj);
                    ngbr(nxt,jj) = 0;
                    nxt = nxtnxt;
                    pthx{cnt}(end+1) = a(nxt,jj+1);
                    pthy{cnt}(end+1) = b(nxt,jj+1);                
                    jj = jj+1;
                else
                    tst = 0;
                end
            end
            cnt = cnt+1;                
        end
    end
end
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