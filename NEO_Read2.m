function [img, head] = NEO_Read2( name, img_num )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if (strcmp(name(end-3:end), '.txt'))
    name = name(1:end-4);
end

%if (nargin > 2)
%    head.Num_Images = num_images;
%end

fid = fopen([name '.txt']);
tmp = fgetl(fid);

if (~strcmp(tmp, '[Dye Settings]'));
    fprintf(1,'\n\n      Warning: File type or version not supported. Aborted.\n');
    img = [];
    head = [];
    return;
end

while(~feof(fid))
    tmp = fgetl(fid);
    if (regexp(tmp, '[\w\s]+ = "?[\w\s\W]+"?'))
        tmp2 = regexp(tmp, ' = ', 'split');
        if (regexp(tmp2{1}(1), '[0-9]'))
            tmp2{1} = ['Val_' tmp2{1}];
        end
        if (regexp(tmp2{2}, '"[\w\s\W]+"'))
            head.(strrep(strrep(tmp2{1},'.',''),' ','_')) = tmp2{2}(2:end-1);
        else
            head.(strrep(strrep(tmp2{1},'.',''),' ','_')) = str2double(tmp2{2});
        end        
    end
end

fclose(fid);
head.AOI_Width = head.ROI_Right - head.ROI_Left + 1;
head.AOI_Height = head.ROI_Top - head.ROI_Bottom + 1;

if (nargin > 1)
    if ((img_num < 1) || (img_num > head.Num_Images))
        fprintf(1,'\n\n      Warning: no image with this number.\n');
        img = [];
    else
        fnum = floor(img_num / 50);
        foffset = img_num - (fnum * 50) - 1;
        pixels = head.AOI_Width * head.AOI_Height;
         
        fid = fopen([name '.' sprintf('%03d', fnum+1)]);
        fseek(fid, pixels * foffset * 2, 'bof');
        img = fread(fid, pixels, 'uint16');
        img = reshape(img, head.AOI_Width, head.AOI_Height)';
    end
else
    img = zeros(head.AOI_Height, head.AOI_Width, head.Num_Images);
    for i = 1:head.Num_Images-1
        [img(:,:,i), head] = PALM_Read( name, i);
    end
end



end

