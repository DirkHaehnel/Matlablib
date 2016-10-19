function err = GridFit(para,y,polyorder)

if nargin<3 | isempty(polyorder)
    polyorder = 1;
end

y = y(:);
x = (1:length(y))'-para(2);

xx = mod(x-1-para(2),para(1))+0.5-para(1)/2;
z = exp(-xx.^2/2/para(3)^2);

for j=1:polyorder
    z(:,1+j) = x.^(j-1);
end

c = z\y;
z = z*c;

plot(x,y,'o',x,z); drawnow

err = sum((y-z).^2);

