%nameFile = 'D:\FASL\SOFI\20110902\Data\';
%DirectoryOut = 'D:\FASL\SOFI\20110902\Tests\';
%actualFile = '06_Q625_GABAR1_X2';

nameFile = strcat(exterNameFile,externActualFile,'.tif');

tic;
tempIm  = imread(nameFile);
dimIm   = size(tempIm);
info = imfinfo(nameFile);
num_images = numel(info);
im = zeros(dimIm(1),dimIm(2),num_images);

averageTimeReadParams = toc;
tic;

im(:,:,1) = tempIm;
clear tempIm dimIm

for j=2:num_images
    im(:,:,j) = imread(nameFile,j, 'Info', info); 
end
timeReadImages = toc;
averageTimeReadImages = timeReadImages/num_images;
clear info num_images j
save(strcat(externDirectoryOut,'Data_',externActualFile,'_OriginalData','.mat'),'-v7.3','im')
%clear im
%return 

tic;
MaximaOriginalImage = max(im,[],3);
averageOriginalImage = mean(im,3);
averageOriginalImage = averageOriginalImage - min(min(averageOriginalImage));
averageOriginalImage = averageOriginalImage./max(max(averageOriginalImage));

MaximaOriginalImage = MaximaOriginalImage - min(min(MaximaOriginalImage));
MaximaOriginalImage = MaximaOriginalImage./max(max(MaximaOriginalImage));

averageTimeOriginalImage = toc;

imwrite(averageOriginalImage,strcat(externDirectoryOut,externActualFile,'_AverageImage.tif'));
imwrite(MaximaOriginalImage,strcat(externDirectoryOut,externActualFile,'_MaxImage.tif'));
tic;
[sof, imAp] = FASLSOFIXAnalysis(im,2);
timeSOFIXAnalysis = toc;
clear im

%return
%[sof2, im2] = SOFIXAnalysis(ch2);
%[sof3, im3] = SOFIXAnalysis(ch3);
%[sof4, im4] = SOFIXAnalysis(ch4);
%[sof12, im12] = SOFIXAnalysis(squeeze(interp1(1:2,permute(cat(4,ch1,ch2),[4 1 2 3]),1.5,'spline')));
%[sof23, im23] = SOFIXAnalysis(squeeze(interp1(1:2,permute(cat(4,ch2,ch3),[4 1 2 3]),1.5,'spline')));
%[sof34, im34] = SOFIXAnalysis(squeeze(interp1(1:2,permute(cat(4,ch3,ch4),[4 1 2 3]),1.5,'spline')));

j=1;
%stck2 = cat(3,sof1(:,:,j),sof12(:,:,j),sof2(:,:,j),sof23(:,:,j),sof3(:,:,j),sof34(:,:,j),sof4(:,:,j));
%stck2 = cat(3,sof(:,:,1),sof(:,:,2));
stck2 = sof(:,:,1);
stck2 = stck2/max(stck2(:));

tic;
% Für das Fourier-apodizieren lauten die Paramter: Pixelgröße 133nm, NA=1.2
% H2O-immersion. -> Brechungsindex in der Probe ~1.334

NA = 1.4; % numerical aperture
n1 = 1.334; % ref. index of sample
n = n1; 
n2 = n1;
d1 = [];
d = 0;
d2 = [];
lambda = 0.625; % emission wavelength in micron
mag = 160; % magnification
pix = mag*0.100/4; % virtual pixel size in micron
sofFourier = FASLSOFIXFourier(stck2,NA,n1,n,n2,d1,d,d2,lambda,mag,pix);
sofFourier = sofFourier/max(sofFourier(:));
clear NA n1 n n2 d1 d d2 lambda mag pix
for j=1:size(sof,3)
    sof(:,:,j)=sof(:,:,j)-min(min(sof(:,:,j)));
    sof(:,:,j)=sof(:,:,j)./max(max(sof(:,:,j)));
    imwrite(sof(:,:,j),strcat(externDirectoryOut,externActualFile,'_Sofi',int2str(j),'.tif'))
    
    sofFourier(:,:,j)=sofFourier(:,:,j)-min(min(sofFourier(:,:,j)));
    sofFourier(:,:,j)=sofFourier(:,:,j)./max(max(sofFourier(:,:,j)));
    imwrite(sofFourier(:,:,j),strcat(externDirectoryOut,externActualFile,'_SofiFourier',int2str(j),'.tif'))
end
save(strcat(externDirectoryOut,'Data_',externActualFile,'_ProccesedData','.mat'))
clear sof sofFourier

return
