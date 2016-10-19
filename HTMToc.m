s = {};
names = dir('*.m');
for j=1:length(names) s{j} = strcat('<a href="',names(j).name,'">',names(j).name,'</a><br>'); end
jj = length(s);
names = dir('*.c');
for j=1:length(names) s{jj+j} = strcat('<a href="',names(j).name,'">',names(j).name,'</a><br>'); end
jj = length(s);
names = dir('*.dll');
for j=1:length(names) s{jj+j} = ['<a href="',names(j).name,'">',names(j).name,'</a><br>']; end

fid = fopen('index.htm','w');

fprintf(fid,'%s','<html>');
fprintf(fid,'%s','<head>');
fprintf(fid,'%s','<title>J&ouml;rg Enderlein''s Matlab Directory</title>');
fprintf(fid,'%s','</head>');
fprintf(fid,'%s','<BODY BACKGROUND="..\back.jpg">');
fprintf(fid,'%s','<h1>');
fprintf(fid,'%s','<center>');
fprintf(fid,'%s','J&ouml;rg Enderlein''s Matlab Directory');
fprintf(fid,'%s','</center>');
fprintf(fid,'%s','</h1>');
fprintf(fid,'%s','<hr>');
fprintf(fid,'%s','<br>');

for j=1:length(s)
    fprintf(fid,'%s',s{j});
end

fclose(fid);