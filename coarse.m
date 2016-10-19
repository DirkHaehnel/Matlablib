function y = Coarse(x,mm);

if prod(size(mm))==1
    mx = mm;
    my = mm;
end

[m,n] = size(x);
m = floor(m/mx)*mx;
n = floor(n/my)*my;

x = x(1:m,1:n);

x = reshape(mean(reshape(x,mx,m*n/mx)),m/mx,n);
y = reshape(mean(reshape(x',my,m*n/mx/my)),n/my,m/my)';


