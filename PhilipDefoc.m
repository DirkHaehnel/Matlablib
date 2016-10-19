% fnames = dir('c:\Joerg\Doc\Microscopy\Defocused\Tinnefeld\*.asc');
% for j=1:100 
%     z(:,:,j)=load(['c:\Joerg\Doc\Microscopy\Defocused\Tinnefeld\' fnames(j).name]); 
% end
% z = z(1:end-1,3:end,:);
% tst = mean(z(:,:,20:end),3);

% z = load('c:\Joerg\Doc\Microscopy\Defocused\Tinnefeld\single strand\3.raw');
% z = load('c:\Joerg\Doc\Microscopy\Defocused\Tinnefeld\double strand\3.raw');

z = load('c:\Joerg\Doc\Microscopy\Defocused\Tinnefeld\22s_1.raw');

ind = ceil(log10(z(:,2:2:end)));ind(isinf(ind))=1;
z = z(:,1:2:end)+z(:,2:2:end)./10.^ind;
im = z(3:end,1:end-2);
[x,y] = meshgrid(1:size(im,2),1:size(im,1));
a = polyfit(x(1,:),mean(im,1),1);
b = polyfit(y(:,1),mean(im,2),1);
im = im-a(1)*x-b(1)*y;

clear model
defocv = (0.8:0.01:1);
for j=1:length(defocv)
    defoc = defocv(j);
    model{j} = PatternGeneration(1.35, 1.51, 1.33, 1.33, [], 0, [], 0.57, 60, defoc, 60*0.12, [], [], [], 10);
end
for j=1:length(defocv)
    [err, bim, cim, sim, xc, yc, bc, cc, sc, len, imm(:,:,j)] = FindPattern(im, model{j}.mask);
end
for j=1:length(defocv)
    tmp=0*model{1}.mask(:,:,1);
    for k=1:size(model{1}.mask,3)
        tmp = tmp + model{j}.mask(:,:,k)*sin(model{j}.theta(k));
    end
    [err, bim, cim, sim, xc, yc, bc, cc, sc, len, imm(:,:,j)] = FindPattern(im, tmp);
    tst(j) = 0;
    for k=1:length(xc)
        tst(j) = tst(j) + err(xc(k),yc(k));
    end
end
for j=1:length(defocv) mim(cat(3,im,imm(:,:,j))); title(defocv(j)); pause; end

return

j = 11;
tmp = 0*model{1}.mask(:,:,1);
for k=1:size(model{1}.mask,3)
    tmp = tmp + model{j}.mask(:,:,k)*sin(model{j}.theta(k));
end
mask = model{j}.mask;

fnames = dir('c:\Joerg\Doc\Microscopy\Defocused\Tinnefeld\*.raw');

for k=1:length(fnames)
    z = load(['c:\Joerg\Doc\Microscopy\Defocused\Tinnefeld\' fnames(k).name]);
    ind = ceil(log10(z(:,2:2:end)));ind(isinf(ind))=1;
    z = z(:,1:2:end)+z(:,2:2:end)./10.^ind;
    im = z(3:end,1:end-2);
    [x,y] = meshgrid(1:size(im,2),1:size(im,1));
    a = polyfit(x(1,:),mean(im,1),1);
    b = polyfit(y(:,1),mean(im,2),1);
    im = im-a(1)*x-b(1)*y;

    [err, bim, cim, sim, xc, yc, bc, cc, sc, len, imm] = FindPattern(im, mask);
    mim(cat(3,im,imm));
    eval(['print -dpng -r300 ' fnames(k).name(1:end-4) 'Fixed'])

    [err, bim, cim, sim, xc, yc, bc, cc, sc, len, imm] = FindPattern(im, tmp);
    mim(cat(3,im,imm));
    eval(['print -dpng -r300 ' fnames(k).name(1:end-4) 'Iso'])
end

return

sig = 50;
mim(cat(3,im,abs(ifft2(exp(-(x-mean(x(:))).^2/2/sig^2-(y-mean(y(:))).^2/2/sig^2).*fftshift(fft2(im))))))
