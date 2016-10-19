function x = CCDRead(name)

fid = fopen(name);
x = fread(fid,inf,'uint8');
fclose(fid);
x = 255-reshape(x,768,576);
