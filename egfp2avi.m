function cx = egfp2avi(name, nn, bin, conv_flag);

% converting image stacks into movies

if nargin==0
    [filename, pathname] = uigetfile('*.tif', 'Interactive mode data:', 0, 0)
    name = [pathname filename];
end
if nargin<2 | isempty(nn)
    nn = 1;
end
if nargin<3 | isempty(bin)
    bin = 0;
end
if nargin<4 | isempty(conv_flag)
    conv_flag = 20;
end
imlen = length(imfinfo(name));
mov = moviein(imlen/nn);
close all
tmp = double(imread(name,1));
for j=2:imlen
    tmp = tmp + double(imread(name,j));
end
figure
set(gcf,'units','normalized','position',[0 0 0.9 0.9])
mim(log(tmp));
[b,a] = ginput(2);
a = round(a); b = round(b);
tmp = tmp(a(1):a(2),b(1):b(2));
if bin>0
    tmp = mconv2(tmp,disk(bin,1/bin^2));
end
if conv_flag>0
    tmp = tmp - mconv2(tmp,disk(conv_flag));            
end
cx = [min(tmp(:)) max(tmp(:))]/imlen*nn;
cx = [0.5 0.8].*cx;
for j=1:nn
    if mod(j-1,nn)==0
        tmp = double(imread(name,j));
    else
        tmp = tmp + double(imread(name,j));
    end
    if mod(j,nn)==0 
        tmp = tmp(a(1):a(2),b(1):b(2));
        if bin>0
            tmp = mconv2(tmp,disk(bin,1/bin^2));
        end
        if conv_flag>0
            tmp = tmp - mconv2(tmp,disk(conv_flag));            
        end
        mim(tmp);
        caxis(cx)
        mov(:,1) = getframe;        
%        tmpold = tmp;
    end
end
for j=nn+1:imlen
    if mod(j-1,nn)==0
        tmp = double(imread(name,j));
    else
        tmp = tmp + double(imread(name,j));
    end
    if mod(j,nn)==0 
        tmp = tmp(a(1):a(2),b(1):b(2));
        if bin>0
            tmp = mconv2(tmp,disk(bin,1/bin^2));
        end
        if conv_flag>0
            tmp = tmp - mconv2(tmp,disk(conv_flag));            
        end
        mim(tmp);
        caxis(cx)
%        mim(tmp-tmpold);
        mov(:,j/nn) = getframe;
%        tmpold = tmp;
    end
end
movie2avi(mov,[name(1:end-4) 'bin' num2str(bin) '.avi'],'compression','none');
close all

