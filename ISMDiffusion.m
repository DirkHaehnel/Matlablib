dname = 'm:\ISM\2012-01\09\Diffusion20nmFlouroSpheres1300Hz1000er\';
%dname = 'm:\ISM\2012-01\06\Diff\Data\';
fnames = dir([dname '*ISM.mat']);
cutstart = 10;

load([dname fnames(1).name])
%len = length(results.SHR)/32^2-cutstart+1;
%timetrace = zeros(len,length(fnames));

len = length(results.SHR)/32^2;
ac = zeros(len,32,32);

for j=1:length(fnames)
    load([dname fnames(j).name])
    im = reshape(results.SHR,32,32,length(results.SHR)/32^2);
    %timetrace(:,j) = squeeze(sum(sum(im(:,:,cutstart:end),1),2))';
    for kx=1:32
        for ky=1:32
            ac(:,kx,ky) = Autocorr(squeeze(im(kx,ky,:))');
        end
    end
    %plot(timetrace(:,j));
    %drawnow
end
%y = Autocorr(timetrace);
close; 
p = Simplex('Rigler',[10 100],[],[],[],[],2:len/2,ac(2:len/2,16,16),0,1)




