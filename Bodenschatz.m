im = tiff2mat('m:\TIRF - Spinning Disc\Dicty 11-12-13\TIRF_Alexa 488_0007.tif');
im = double(im);
mim(im(:,:,1));
[a,b] = PickVec;

for j=1:20 
    [err, bim, cim, sim, xc, yc, bc, cc, sc, len, imm] = FindPattern(im(a,b,1),Disk(3,['exp(-rad.^2/(1+0.1*' num2str(j) '))'])); 
    mim(cat(3,im(a,b,1),imm)); title(mint2str(j,2)); 
    pause; 
end

imm = im(a,b,:);
for j=1:size(im,3)
    [err, bim, cim, sim, xc, yc, bc, cc, sc, len, tst] = FindPattern(im(a,b,j),Disk(3,['exp(-rad.^2/2)']));
    if ~isempty(tst)
        imm(:,:,j) = tst;
    else
        imm(:,:,j) = 0;
    end
    mim(cat(3,im(a,b,j),imm(:,:,j))); title(mint2str(j,2)); 
    eval(['print -dpng -r200 tmp' mint2str(j,3)]);
end
    