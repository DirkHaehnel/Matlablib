% Localization and tracking of single molecule images on a surface

clear all
close all
name = 'D:\Joerg\Doc\Patra\EGFP\Cy5step_10000times_water_5s_100seq_a1.tif';
th = 50;
rr = 20;

imlen = length(imfinfo(name));
x=double(imread(name,1)); 
mm = [min(x(:)) max(x(:))];
for j=2:imlen 
    tmp = double(imread(name,j));
    mm(1) = min([mm(1) min(tmp(:))]);
    mm(2) = max([mm(2) max(tmp(:))]); 
    x = x+tmp;
end

imlen = length(imfinfo(name));
nn = 2;
mov = moviein(imlen/nn);
for j=1:nn
    if mod(j-1,nn)==0
        tmp = double(imread(name,j));
    else
        tmp = tmp + double(imread(name,j));
    end
    if mod(j,nn)==0 
        tmp1 = mconv2(tmp,disk(10));            
        mim(tmp-tmp1);
        cx = [0.5 0.8].*[min(min(tmp-tmp1)) max(max(tmp-tmp1))];
        caxis(cx)
        mov(:,j/nn) = getframe;
    end
end
for j=nn+1:imlen
    if mod(j-1,nn)==0
        tmp = double(imread(name,j));
    else
        tmp = tmp + double(imread(name,j));
    end
    if mod(j,nn)==0 
        tmp1 = mconv2(tmp,disk(10));            
        mim(tmp-tmp1);
        caxis(cx)
        mov(:,j/nn) = getframe;
    end
end

movie2avi(mov,[name(1:end-4) 'bin' num2str(nn) '.avi'],'compression','none');

x1 = mconv2(x,disk(1));
x2 = mconv2(x,disk(10));

%[cc,cl] = cluster(max2(mconv2(x,disk(10,0.3))).*(x1>1.01*x2));
[cc,cl] = cluster(x1>1.06*x2);

% if length(cl)>240
%     j1 = length(cl)-240+1;
% else
%     j1=1;
% end
j1 = 1;
h = waitbar(0,'Finding coordinates');
[xx,yy]=meshgrid(1:size(x,2),1:size(x,1));
for j=j1:length(cl) 
    a(j-j1+1)=round(mean(xx(cc==j))); 
    b(j-j1+1)=round(mean(yy(cc==j))); 
    waitbar(j/(length(cl)-j1+1));
end
close(h);
mim(x-x2); hold; plot(a,b,'oy'); hold; drawnow
 
b(a<rr+1) = [];
a(a<rr+1) = [];
a(b<rr+1) = [];
b(b<rr+1) = [];
b(a>size(x,2)-rr)=[];
a(a>size(x,2)-rr)=[];
a(b>size(x,1)-rr)=[];
b(b>size(x,1)-rr)=[];

clear z

h = waitbar(0,'Processing');
for j=1:imlen 
    x = double(imread(name,j)); 
%     x = mconv2(x,disk(1,1));
%     x = x - mconv2(x,disk(10));
    for k=1:length(a)
        z(j,k) = sum(sum(x(b(k)-rr:b(k)+rr,a(k)-rr:a(k)+rr)));
    end
    waitbar(j/imlen,h);
end
close(h);

% colnum=ceil(sqrt(size(z,4))*4/3);
% rownum = ceil(size(z,4)/colnum);
colnum = 10;
rownum = ceil(size(z,4)/colnum);
mov=moviein(size(z,3)); 
cx = [min(z(:)) max(z(:))]*4/5;
for j=1:size(z,3)
    tmp=ones(rownum*(2*rr+1)+rownum-1,colnum*(2*rr+1)+colnum-1)*inf; 
    for k=1:size(z,4)
        tmp((ceil(k/colnum)-1)*(2*rr+2)+1:ceil(k/colnum)*(2*rr+2)-1, ...
            rem(k-1,colnum)*(2*rr+2)+1:(rem(k-1,colnum)+1)*(2*rr+2)-1) = z(:,:,j,k); 
    end; 
    mim(tmp); caxis(cx);
    for k=1:colnum
        text((k-0.5)*(2*rr+2),-rr,int2str(k),'horizontalalignment','center');
    end
    for k=1:rownum
        text(-rr,(k-0.5)*(2*rr+2),-rr,int2str((k-1)*colnum),'horizontalalignment','center');
    end
    mov(:,j)=getframe; 
end

% auto = [];
% for j=1:size(z,3)-1
%     for k=1:size(z,4)
%         auto(j,k) = sum(sum(z(:,:,j,k).*z(:,:,j+1,k)))/sqrt(sum(sum(z(:,:,j+1,k).*z(:,:,j+1,k)))*sum(sum(z(:,:,j,k).*z(:,:,j,k)))); 
%     end; 
% end

break

for k=1:size(z,4)
    cx = [min(min(min(z(:,:,:,k)))) max(max(max(z(:,:,:,k))))]*4/5;
    colnum=ceil(sqrt(size(z,3))*4/3);
    rownum = ceil(size(z,3)/colnum);
    tmp=ones(rownum*(2*rr+1)+rownum-1,colnum*(2*rr+1)+colnum-1)*inf; 
    for j=1:size(z,3)
        tmp((ceil(j/colnum)-1)*(2*rr+2)+1:ceil(j/colnum)*(2*rr+2)-1, ...
            rem(j-1,colnum)*(2*rr+2)+1:(rem(j-1,colnum)+1)*(2*rr+2)-1) = z(:,:,j,k); 
    end
    mim(tmp); caxis(cx); times16; title(['image # ' int2str(k)]); 
    pause
end

