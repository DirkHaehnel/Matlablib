cd d:/joerg/doc/patra/sauerdots

close all
clear all

load picks
mm = (size(imm,1)-1)/2;
ind = 1:mm;
back = disk(mm);
pos = 1:36;
for k=1:size(imm,3)
    im = imm(:,:,k);
    for j=1:36
        tmp = back.*imrotate(im,(j-1)*5,'bicubic','crop');
        joerg(1,j) = ind*sum(tmp(ind,:).^4')'/sum(sum(tmp(ind,:).^4));
        joerg(2,j) = ind*sum(tmp(size(tmp,1)-ind,:).^4')'/sum(sum(tmp(size(tmp,1)-ind,:).^4));    
    end
    j = pos(sum(joerg)==max(sum(joerg)))
    imr(:,:,k) = imrotate(im,(j-1)*5,'bicubic','crop');
    subplot(7,7,k); mim(imr(:,:,k)); drawnow
end

