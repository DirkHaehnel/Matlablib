% Localization and tracking of single molecule images on a surface

clear all
close all
name = 'd:/joerg/doc/patra/egfp/eGFP_filterstock_17mW_5s_seq60_a1.tif';
th = 50;
rr = 20;

imlen = length(imfinfo(name));
x=double(imread(name,1)); 
for j=2:imlen 
    x=x+double(imread(name,j)); 
end

x1 = mconv2(x,disk(1));
x2 = mconv2(x,disk(10));

[cc,cl] = cluster(max2(mconv2(x,disk(10,0.3))).*(x1>1.01*x2));

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
% 
% r=sqrt((ones(length(a),1)*a-a'*ones(1,length(a))).^2+(ones(length(b),1)*b-b'*ones(1,length(b))).^2);
% j=1; k=2;
% while j<length(a)
%     while k<=length(a)
%         if r(j,k)<rr
%             a(j) = round(mean(a([j k])));
%             b(j) = round(mean(b([j k]))); 
%             a(k) = [];
%             b(k) = [];
%             r(:,k) = [];
%             r(k,:) = [];
%             r(:,j) = sqrt((a(j)-a).^2+(b(j)-b).^2)';
%             r(j,:) = r(:,j)';
%         else
%             k = k+1;
%         end
%     end
%     j = j+1;
%     k = j+1;
% end
%             
b(a<rr+1) = [];
a(a<rr+1) = [];
a(b<rr+1) = [];
b(b<rr+1) = [];
b(a>size(x,2)-rr)=[];
a(a>size(x,2)-rr)=[];
a(b>size(x,1)-rr)=[];
b(b>size(x,1)-rr)=[];
% clear z

h = waitbar(0,'Processing movie');
for j=1:imlen 
    x = double(imread(name,j)); 
    x = mconv2(x,disk(1,1));
    x = x - mconv2(x,disk(10));
    for k=1:length(a)
        z(:,:,j,k) = x(b(k)-rr:b(k)+rr,a(k)-rr:a(k)+rr);
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

