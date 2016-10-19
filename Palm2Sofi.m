function [sofi_im, mean_im, im4] = Palm2Sofi(palm_im)

mean_im = mean(palm_im,3);
palm_im = palm_im - repmat(mean_im,[1 1 size(palm_im,3)]);
for m = 1:3
    tmp = palm_im(:,:,1:end-m).*palm_im(:,:,1+m:end);
    sofi_im(:,:,m,1) = mean(tmp,3); 
    tmp = tmp(:,:,1:end-m).*palm_im(:,:,1+2*m:end);
    sofi_im(:,:,m,2) = mean(tmp,3); 
    tmp = tmp(:,:,1:end-m).*palm_im(:,:,1+3*m:end);
    sofi_im(:,:,m,3) = mean(tmp,3); 

%     -10 C[2] C[3] + DR[5]
%     -15 C[2]^3 - 10 C[3]^2 - 15 C[2] C[4] + DR[6]
%     -105 C[2]^2 C[3] - 35 C[3] C[4] - 21 C[2] C[5] + DR[7]
%     -105 C[2]^4 - 210 C[2]^2 C[4] - 35 C[4]^2 - 56 C[3] C[5] - 28 C[2] (10 C[3]^2 + C[6]) + DR[8]    
end

im4 = abs(sofi_im(:,:,1,3)-sofi_im(:,:,1,1).^2-sofi_im(:,:,2,1).^2-sofi_im(:,:,3,1).^2);
close; mim(cat(4,cat(3,mean_im,sofi_im(:,:,1,1)),cat(3,abs(sofi_im(:,:,1,2)),im4)));

