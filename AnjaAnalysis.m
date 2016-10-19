clear all
close all

dirname = 'C:\Users\Anja Huss\Pictures\Daten\20012012 vergleich\';
fnames = dir([dirname '*.tif']);
nn = length(fnames);

if 1 % convert tiff data to mat-files
    for j=1:nn
        info =imfinfo([dirname fnames(j).name]);
        len = numel(info);
        for k=1:len
            im(:,:,k) = imread([dirname fnames(j).name], k, 'Info', info);
            disp(k)
        end
        eval(['save ''' dirname fnames(j).name(1:end-4) ''' im'])
    end
end


if 1 % check decay times
    for j=1:nn
        eval(['load ''' dirname fnames(j).name(1:end-4) '''']);
        im = double(im);
        ac = squeeze(sum(sum(abs(fft(im-repmat(mean(im,3),[1 1 size(im,3)]),[],3)).^2,1),2));
        close all
        p(:,j) = Simplex('ExpFun',[1 2],[],[],[],[],1:20,ac(2:21),2);
        [err, c(:,j)] = ExpFun(p(:,j), 1:20, ac(2:21));
        tau_mean(j) = (c(2:3,j)'*p(:,j))/sum(c(2:3,j));
    end
end

