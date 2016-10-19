

dirname = 'Z:\SOFI\Neurons Data\2012-12 samples\2012-12-30_1 Neurons\2012-12-30_1 raw SOFI\Test FR\';

fnames = dir([dirname '*9.mat']);
    for j=1:length(fnames)
        load([dirname fnames(j).name]);
        res = abs(fft2(fftshift(ifftshift(ifft2(final)).*weight.*filter)));   % inside the most inner bracket is the name of the variable save in the mat file
        save([dirname fnames(j).name '_res.txt'],'res','-ascii');
    end
    