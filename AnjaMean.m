clear all 
close all

[filename, pathname] = uigetfile('*.mat', 'Pick an mat-File','\\jesrv\WG-Data\*');
if isequal(filename,0)
    return
else
    
   load(fullfile(pathname, filename));
   
  

end
mim(mean(img,3));

%wenn fehler schauen was im workspace steht evtl img in im oder data
%wandeln.