function [pos, z] = Myosin(filename,av,bv)

if nargin==0
    [filename, pathname] = uigetfile('*.tif;*.tiff', 'Interactive mode data:', 0, 0)
    name = {[pathname filename]}; 
else
    if iscell(filename)
        name = filename;
    else
        name = {filename};
    end
end

if nargin<2
    for filecnt=1:length(name)
        tmp = mtiffread(name{filecnt});
        clear im
        for j=1:length(tmp) im(:,:,j)=double(tmp(j).data); end
        clf; mim(sum(im,3));
        [b,a] = ginput(2);
        a=round(a); b=round(b);
        a(1) = max([a(1) 1]);
        a(end) = min([a(end) size(im,1)]);
        b(1) = max([b(1) 1]);
        b(end) = min([b(end) size(im,2)]);
        av(:,filecnt) = a; bv(:,filecnt) = b;
    end
end

[xx,yy] = meshgrid(-2:2,-2:2);
[xi,yi] = meshgrid(-2:0.1:2,-2:0.1:2);
for filecnt=1:length(name)
    tmp = mtiffread(name{filecnt});
    outfile = findstr(name{filecnt},'.');
    outfile = name{filecnt}(1:outfile(end)-1);
    for j=1:length(tmp) im(:,:,j)=double(tmp(j).data); end
    a = av(:,filecnt); b = bv(:,filecnt);
    close all
    for j=1:size(im,3)
        z(:,:,j) = im(a(1):a(2),b(1):b(2),j);
        [pos(1,j) pos(2,j) pos(3,j), pos(4,j)] = GaussEllipse([],[],z(:,:,j));
    end
end

