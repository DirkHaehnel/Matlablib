n = 1/1.5;
x = tan(pi/8);

q = (1-n^2)/2/x;
w1 = q - sqrt(q^2-n^2);
w2 = q + sqrt(q^2-n^2);

phi1 = acos(sqrt((1-n^2)/(1+w1^2)));
phi2 = acos(sqrt((1-n^2)/(1+w2^2)));

theta = asin(n):pi/1e3:pi/2;
w = sqrt((1-n^2)./cos(theta).^2-1); 
plot(theta/pi*180,2*(atan(w/n^2)-atan(w))/pi*180,theta/pi*180,ones(size(theta))*45,':')

