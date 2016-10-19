clear all 
close all

[filename, pathname] = uigetfile('*.mat', 'Pick an ISM-File','\\jesrv\WG-Data\ISM\*');
if isequal(filename,0)
    return
else
    load(fullfile(pathname, filename));
   
  

end

mim(sum(img,3));