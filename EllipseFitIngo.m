function err = EllipseFit(para,im)

% para(1) = x-Koordinate Mittelpunkt
% para(2) = y-Koordinate Mittelpunkt
% para(3) = Exzentrizität
% para(4) = Winkel der großen Achse zur x-Achse
% para(5) = min. Radius 
% para(6) = max. Radius

im = zeros([100 100]);

para = [50 50 0.7 pi/4 15 20];

[m,n] = size(im);

ez    = sqrt(1-para(3)^2);
theta = para(4);

[x,y] = meshgrid((1:m)-para(1),para(2)-(1:n));

ra = sqrt(x.^2+y.^2);

phi = zeros(size(im));
phi((x>=0)&(y>=0)) =        atan((y((x>=0)&(y>=0)))./(x((x>=0)&(y>=0))));
phi((x< 0)&(y>=0)) =   pi + atan((y((x< 0)&(y>=0)))./(x((x< 0)&(y>=0))));
phi((x< 0)&(y< 0)) =   pi + atan((y((x< 0)&(y< 0)))./(x((x< 0)&(y< 0))));
phi((x>=0)&(y< 0)) = 2*pi + atan((y((x>=0)&(y< 0)))./(x((x>=0)&(y< 0))));

% xx  = (sin(phi).*cos(theta) - ez.*cos(phi).*sin(theta));
% yy  = (sin(phi).*sin(theta) + ez.*cos(phi).*cos(theta));

xx  = ( sin(phi).*cos(theta) - ez.*cos(phi).*sin(theta) );
yy  = ( sin(phi).*sin(theta) + ez.*cos(phi).*cos(theta) );

zz = sqrt((xx.^2+yy.^2));

re = ra./zz;
%re(re<para(5)) = 0;
%re(re>para(6)) = 0;
mim(re);

% z = z>=para(3) & z<=para(4);
% col = ones(prod(size(im)),1);
% c = [col z(:)]\im(:);
% 
% mim(z.*im); drawnow;
% 
% err = sum((im(:) - [col z(:)]*c).^2);
