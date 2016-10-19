function x = alvread(name);

fin = fopen(name);
for j=1:17 
    fgetl(fin);
end

x = zeros(183,12);
for j=1:183
    tmp = fgetl(fin);
    x(j,:) = str2num(tmp);
end