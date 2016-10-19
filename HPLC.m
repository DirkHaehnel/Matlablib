function [pos, area] = HPLC(name)

if nargin==0
    [filename, pathname] = uigetfile('*.arw')
    name = [pathname filename];
end 

% name = 'D:\Daten\NovakJens\Jens Novak GC-Assay43629.arw'

fin = fopen(name,'r');

for j=1:2 fgetl(fin); end

cnt = 1;
while ~feof(fin)
    tmp = fgetl(fin);
    ind = findstr(tmp,',');
    tmp(ind) = '.';
    x(cnt,:) = str2num(tmp);
    cnt = cnt+1;
end
fclose(fin);
ind = x(:,1)>6 & x(:,1)<8;
x = x(ind,:);
x(:,2) = x(:,2)-min(x(:,2));
plot(x(:,1),x(:,2))
pos = ginput;
pos = pos(:,1);
ind = x(:,1)>min(pos) & x(:,1)<max(pos);
p = simplex('Lorenz',[pos(2:end-1); 0.1*ones(size(pos,1)-2,1)], [min(pos)*ones(size(pos,1)-2,1); zeros(size(pos,1)-2,1)], [max(pos)*ones(size(pos,1)-2,1); inf*ones(size(pos,1)-2,1)], [], [], x(ind,1), x(ind,2), 3);

for j=1:5
    p = simplex('Lorenz',p, [min(pos)*ones(size(pos,1)-2,1); zeros(size(pos,1)-2,1)], [max(pos)*ones(size(pos,1)-2,1); inf*ones(size(pos,1)-2,1)], [], [], x(ind,1), x(ind,2), 3);
end

[err, c, z, A] = Lorenz(p,x(ind,1), x(ind,2));

plot(x(:,1),x(:,2),x(ind,1),A.*(ones(size(A,1),1)*c'))

pos = p(1:end/2);
area = pi*p(end/2+1:end).*c(1:length(pos));