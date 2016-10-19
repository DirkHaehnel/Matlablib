function DeInterlace(name)

file_list = dir([name '*.tif']);
len = length(name)+1;
j = 1;
while j<=length(file_list)
    if strncmp([name '.'], file_list(j).name, len) file_list(j)=[]; end
    j = j+1;
end
    
nmax = floor(log10(2*length(file_list))) + 1;
for j=1:length(file_list)
    im = imread(file_list(j).name);
    ind = repmat(1:2:size(im,1),2,1);
    ind = reshape(ind,prod(size(ind)),1);    
    imwrite(im(ind,:),[name 'NI' mint2str(2*j-1,nmax) '.tif'],'tif', 'compression', 'none');
    ind = repmat(2:2:size(im,1),2,1);
    ind = reshape(ind,prod(size(ind)),1);    
    imwrite(im(ind,:),[name 'NI' mint2str(2*j,nmax) '.tif'],'tif', 'compression', 'none');    
end

