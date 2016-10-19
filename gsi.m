function [ims, imms] = gsi(im,imm);

% function gsi picks a subimage from a displayed one

mim(im);
[y,x] = ginput(2);
x = round(x); 
y = round(y);

ims = im(x(1):x(2),y(1):y(2));
if nargin>1
    imms = imm(x(1):x(2),y(1):y(2));
end
