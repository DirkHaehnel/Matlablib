% program for the calculation of the detection efficiency distribution
% infinitely long slit with width d

NA = 0.85;
psi = asin(NA/1.33);
d = 10;
k1 = 1/sin(psi);
k2 = k1*cos(psi);

close all
y = 0:0.2:10;
clear res;

for x = 0:40
	if x==0
		tst1 = -psi*(-d/2<y).*(y<d/2);
		tst2 = psi*(-d/2<y).*(y<d/2);
	else
		tst1 = max([atan(-(y + d/2)/abs(x)); -psi*ones(size(y))]);
		tst2 = min([atan((d/2 - y)/abs(x)); psi*ones(size(y))]);
	end
	delta1 = sqrt(1 - k1^2*sin(tst1).^2) + 1e-10;
	delta2 = sqrt(1 - k1^2*sin(tst2).^2) + 1e-10;
	res = [res; k1*asin(k1*sin(tst2)) - k2*atan(k2*sin(tst2)./delta2) - k1*asin(k1*sin(tst1)) + k2*atan(k2*sin(tst1)./delta1)];
end
res = res.*(imag(res)==0);
res = 2*sin(psi)/(2*pi*(1-cos(psi)))*res;
if 1
	tmp = res;
	tmp(:,1) = [];
	res = [fliplr(tmp) res];
	tmp = res;
	tmp(1,:) = [];
	res = [flipud(tmp); res];
end;
% close; surf(-10:0.2:10, -40:40,res); view([-37.5,50]); colormap('jet');
title('Collection Efficiency Function', 'FontName', 'Times', 'FontSize', 24);
close; res(41,26) = 0.5; mesh(-10:0.2:10, -40:40,res); view([-40,50]); colormap([zeros(64,1), zeros(64,1), zeros(64,1)])
% colormap('gray'); brighten(0.9);
xlabel('y in micrometer');
ylabel('x in micrometer');
zlabel('CEF');
times16;

