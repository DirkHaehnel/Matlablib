% Localization and tracking of single molecule images on a surface

clear all
close all
cd 'D:\Joerg\Doc\Patra\Dendrimer\'

name=  'G1R4_5percent_1s_seq100_31mW_df1m.tif'; 
%name = 'G1R4_5percent_1s_seq300_31mW_df1m.tif';
imlen = length(imfinfo(name));
jlen = 8;
klen = 8;
xx = double(imread(name)); 
[m,n] = size(xx);
for j=1:jlen
    for k=1:klen
        x = xx((1:m/jlen)+(j-1)*m/jlen, (1:n/klen)+(k-1)*n/klen);
%       x = mconv2(x,disk(1,1));
        x = x - mconv2(x,disk(10));
        cx(j,k,:) = 0.8*[min(x(:)) max(x(:))];
    end
end
h=figure;
%set(h,'position',[10 40 1000 640],'DoubleBuffer','on');
set(h,'position',[10 40 3*n/klen 3*m/jlen],'DoubleBuffer','on');
for j=1:jlen
    for k=1:klen
        x = xx((1:m/jlen)+(j-1)*m/jlen, (1:n/klen)+(k-1)*n/klen);
%       x = mconv2(x,disk(1,1));
        x = x - mconv2(x,disk(10));
        aviobj{j,k} = avifile([name(1:end-4) mint2str(j,1) mint2str(k,1) '.avi'],'Compression','none');        
        pcolor(x); axis off; caxis(squeeze(cx(j,k,:))); 
        axis image; shading flat; colormap hot; drawnow
        frame = getframe(gca);
        aviobj{j,k} = addframe(aviobj{j,k},frame);
    end
end
for r=2:imlen 
    xx = double(imread(name,r)); 
    for j=1:jlen
        for k=1:klen
            x = xx((1:m/jlen)+(j-1)*m/jlen, (1:n/klen)+(k-1)*n/klen);
%            x = mconv2(x,disk(1,1));
            x = x - mconv2(x,disk(10));
            pcolor(x); axis off; caxis(squeeze(cx(j,k,:))); 
            axis image; shading flat; colormap hot; drawnow
            frame = getframe(gca);
            aviobj{j,k} = addframe(aviobj{j,k},frame);
        end
    end
end
for j=1:jlen
    for k=1:klen
        aviobj{j,k} = close(aviobj{j,k});
    end
end
