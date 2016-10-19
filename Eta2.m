% program for the calculation of the efficient detection volume
% rectangular slit d*l

clear
NA = 0.55;
psi = asin(NA/1.33);
tpsi = tan(psi);
tpsi2 = tan(psi)^2;
d = 10;
l = 10;
w1 = 25; % elliptical x-axis of laser beam
w2 = 5; % elliptical z-axis of laser beam

x0 = 0:50;
y0 =0:100;
z0 = 0:50;
fac = 1/pi/(1-cos(psi))/2;
int1 = exp(-2*x0.*x0/w1^2)'*ones(size(y0));
int2 = exp(-2*x0.*x0/w2^2)'*ones(size(y0));
vol1 = 0;
vol2 = 0;

for k= 1:length(z0)
	z = z0(k);
for i = 1:length(x0)
		x = x0(i);
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
	vol1 = vol1 + sum(sum(res.*int1))*exp(-2*z*z/w2^2)
	vol2 = vol2 + sum(sum(res.*int2))*exp(-2*z*z/w1^2)
	vol1/vol2
end
