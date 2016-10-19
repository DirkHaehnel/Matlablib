% program for the calculation of the detection efficiency distribution
% complete rectangular slit d*l

NA = 0.85;
psi = asin(NA/1.33);
tpsi = tan(psi);
tpsi2 = tan(psi)^2;
d = 10;
l = 10;

close
z = 0;
res = [];
x0 = 0:40;
y0 =0:0.2:10;
fac = 1/pi/(1-cos(psi))/2;

for i = 1:length(x0)
	x = abs(x0(i));
	for j = 1:length(y0)
		y = y0(j);
		if x==0
			res(i,j) = (-d/2<y)*(y<d/2)*(-l/2<z)*(z<l/2);
		else
			t1 = max(-d/2, y-abs(x)*tpsi);
			t2 = min(d/2, y+abs(x)*tpsi);
			step = (t2-t1)/100;
			s = t1 + ((1:100)-0.5)*step;
			f1 = max([-l/2*ones(size(s)) - z; - sqrt(x^2*tpsi2-(s-y).^2)]);
			f2 = min([l/2*ones(size(s)) - z; sqrt(x^2*tpsi2-(s-y).^2)]);
			res(i,j) = fac*step*sum(abs(x)./(x^2+(s-y).^2).*(f2./sqrt(x^2 + f2.^2 + (s-y).^2) - f1./sqrt(x^2 + f1.^2 + (s-y).^2)));
		end
	end
end
res = res.*(imag(res)==0);
res = res.*(res>0);
if 1
	tmp = res;
	tmp(:,1) = [];
	res = [fliplr(tmp) res];
	tmp = res;
	tmp(1,:) = [];
	res = [flipud(tmp); res];
end
close; surf(-10:0.2:10, -40:40,res); view([-37.5,50]); colormap('gray'); brighten(0.9);
set(gca,'FontName','Times');
xlabel('y in micrometer','FontName','Times','Rotation',20);
ylabel('x in micrometer','FontName','Times','Rotation',-35);
zlabel('CEF','FontName','Times');

