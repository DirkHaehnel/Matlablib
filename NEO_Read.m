function [img, head] = NEO_Read( directory, width, height, img_num )
% Routine for reading Neo sCMOS spooled .dat files
% directory : directory where dat files are
% width : image width in pixels (might need to add some blank pixels)
% hight : image hight in pixels
% img_num : num of image to read

name = [directory '\acquisitionmetadata.ini'];
if (~exist(name , 'file') || ~exist([directory '\0000000000spool.dat'], 'file'))
    fprintf(1,'\n\n      Warning: File not found.\n');
    img = [];
    head = [];
    return;
end

if (nargin < 3)
    head.AOI_Height = 750;
    head.AOI_Width = 750;
else
    head.AOI_Height = height;
    head.AOI_Width = width;
end

fid = fopen(name);
tmp = fgetl(fid);

if (~strfind(tmp, '[data]'));
    fprintf(1,'\n\n      Warning: File type or version not supported. Aborted.\n');
    img = [];
    head = [];
    return;
end

fnames = dir([directory '\*spool.dat']);

while(~feof(fid))
    tmp = fgetl(fid);
    if (regexp(tmp, '[\w\s]+ = "?[\w\s\W]+"?'))
        tmp2 = regexp(tmp, ' = ', 'split');
        if (regexp(tmp2{1}(1), '[0-9]'))
            tmp2{1} = ['Val_' tmp2{1}];
        end
        if (regexp(tmp2{2}, '"[\w\s\W]+"'))
            head.(strrep(tmp2{1},' ','_')) = tmp2{2}(2:end-1);
        else
            head.(strrep(tmp2{1},' ','_')) = str2double(tmp2{2});
        end        
    end
end

head.Num_Images = length(fnames) * head.ImagesPerFile;

fclose(fid);

if (nargin > 3)
    if ((img_num < 1) || (img_num > head.Num_Images))
        fprintf(1,'\n\n      Warning: no image with this number.\n');
    else
        fnum = floor(img_num / head.ImagesPerFile) + 1;
        foffset = img_num - ((fnum - 1) * head.ImagesPerFile) - 1;
        pixels = head.AOI_Width * head.AOI_Height;
         
        fid = fopen([directory '\' fnames(fnum).name ]);
        fseek(fid, head.ImageSize * foffset * 2, 'bof');
        img = fread(fid, pixels, 'uint16', 0 , 'l');
        img = reshape(img, head.AOI_Width, head.AOI_Height)';
    end
else
    img = zeros(head.AOI_Height, head.AOI_Width, head.Num_Images);
    for i = 1:head.Num_Images-1
        [img(:,:,i), head] = PALM_Read( name, i);
    end
end



end

