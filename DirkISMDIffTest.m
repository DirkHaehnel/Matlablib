% program for OTF determination

clear all 
close all


[filename, pathname] = uigetfile('*.mat', 'Pick an ISM-File','\\jesrv\WG-Data\ISM\*');
if isequal(filename,0)
    return
else
    load(fullfile(pathname, filename));
end


imraw=reshape(results.SHR,32,32,size(results.SHR));

mim(squeeze((sum(imraw,3))));

