clear all
close all

dirname = 'C:\Users\Anja Huss\Pictures\Daten\20012012 vergleich\';
fnames = dir([dirname '*.tif']);

nn = length(fnames);

% tst = zeros(1,nn);
% for j=1:nn
%     for k=2:5
%         if strfind(fnames(j).name,['_X' int2str(k) '.tif']);
%             tst(j) = 1;
%         end
%     end
% end
% fnames(logical(tst)) = [];

if 1 % convert tiff data to mat-files
    for j=1:nn
        info =imfinfo([dirname fnames(j).name]);
        len = numel(info);
        for k=1:len
            im(:,:,k) = imread([dirname fnames(j).name], k, 'Info', info);
            disp(k)
        end
        %         for jj=2:5
        %             if exist([dirname fnames(j).name(1:end-4) '_X' int2str(jj) '.tif'])
        %                 kmax = k;
        %                 info =imfinfo([dirname fnames(j).name(1:end-4) '_X' int2str(jj) '.tif']);
        %                 len = numel(info);
        %                 for k=kmax+1:kmax+len
        %                     im(:,:,k) = imread([dirname fnames(j).name(1:end-4) '_X' int2str(jj) '.tif'], k-kmax, 'Info', info);
        %                     disp(k)
        %                 end
        %             end
        %         end
        eval(['save ''' dirname fnames(j).name(1:end-4) '''.mat im'])
    end
end

if 0 % check decay times
    for j=1:nn
        eval(['load ''' dirname fnames(j).name(1:end-4) '''' '.mat']);
        im = double(im);
        mm = repmat(mean(im,3),[1 1 100])
        ac = squeeze(sum(sum(abs(fft(im(:,:,1:100)-mm,[],3)).^2,1),2));
        close all
        p(:,j) = Simplex('ExpFun',[1 2],[],[],[],[],1:20,ac(2:21),2);
        [err, c(:,j)] = ExpFun(p(:,j), 1:20, ac(2:21));
        tau_mean(j) = (c(2:3,j)'*p(:,j))/sum(c(2:3,j));
    end
end

if 1 % check decay times (cross-correlation)
    p=zeros(2,nn,10);
    tau_mean=zeros(nn,10);
    for j=1:nn
        cc=0;
        for k=1:10
            eval(['load ''' dirname fnames(j).name(1:end-4) '''' '.mat']);
            im = double(im);
%             im=im(1:40,1:40,:);
            mm = repmat(mean(im(:,:,(k-1)*100+(1:100)),3),[1 1 100]);
            cc = cc+0.1*squeeze(sum(sum(real(fft(im(1:end-1,:,(k-1)*100+(1:100))-mm(1:end-1,:,:),[],3).*conj(fft(im(2:end,:,(k-1)*100+(1:100))-mm(2:end,:,:),[],3))),1),2));
            close all
            p(:,j,k) = Simplex('ExpFun',[1 2],[],[],[],[],0:19,cc(1:20),2);
            [err, c(:,j)] = ExpFun(p(:,j), 0:19, cc(1:20));
            tau_mean(j,k) = (c(2:3,j)'*p(:,j))/sum(c(2:3,j));
        end
    end
end

