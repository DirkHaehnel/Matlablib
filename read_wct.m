function [head, data] = read_wct(name);

close all

fid = fopen(name,'r');

if fid ~= -1
   t = '';
   
   while strcmp(t,'=')==0
     t = fgets(fid,1);
   end;
   t = fgets(fid,5);
   width = str2num(t);
   
   while strcmp(t,'=')==0
     t = fgets(fid,1);
   end;
   t = fgets(fid,5);
   height = str2num(t);
   
   while strcmp(t,'=')==0
     t = fgets(fid,1);
   end;
   t = fgets(fid,5);
   xsize = str2num(t);
   
   while strcmp(t,'=')==0
     t = fgets(fid,1);
   end;
   t = fgets(fid,5);
   ysize = str2num(t);

   while strcmp(t,';')==0
     t = fgets(fid,1);
   end;
   t = fgets(fid,2);
   
   data = [];  
   for i = 1:height
     tmp = [];
     [tmp, count ] = fscanf(fid,'%d%c',[2 width]); 
     t = fgets(fid,3);
     data = [data; tmp(1,:)];   

   end;   
   fclose(fid);
   
   head = struct('widhth', width,  ...
                 'height', height, ...
                 'x_size', xsize,  ...
                 'y_size', ysize);
end
