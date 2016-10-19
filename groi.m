function [ims, imms] = groi(im,imm)

% function groi picks a subimage from a displayed one

mim(im);
if nargin>1
    if length(imm)==1
        [y,x] = ginput(1);
        x = round(x); 
        y = round(y);
        ims = im(x-imm:x+mm,y-mm:y+mm);
    else
        [y,x] = ginput(2);
        x = round(x); 
        y = round(y);
        ims = im(x(1):x(2),y(1):y(2));        
        imms = imm(x(1):x(2),y(1):y(2));
    end
else
    [y,x] = ginput(2);
    x = round(x); 
    y = round(y);
    ims = im(x(1):x(2),y(1):y(2));        
end
