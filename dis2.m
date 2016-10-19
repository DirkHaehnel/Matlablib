% Programm zur Bifurkationsuntersuchung eines katalytischen Festbettes

a = 0.5;
b = 0.5;
c = 7;

len = 20
len2 = 2*len+1;
if ~(exist('radii')==1)
	radii = zeros(3, len2);
	k3 = -len:len;
	count = 0;
	for k1 = -len:len
	k1
		for k2 = -len:len
			radii(:, count*len2+1:(count+1)*len2) = [2*k1+k2+k3; sqrt(3)*k2+k3/sqrt(3); sqrt(8/3)*k3];
			count = count + 1;
		end
	end
	radii = sqrt(sum(radii.*radii));
	radii(radii==0)=[];
end

% save radii.mat radii
close
clear null
for jj=1:5
	r = 3+ jj;
	for ii=1:100
		phi=ii*pi/100;
		omega = sqrt(1+i*phi);
		c0 = omega*r*cosh(omega*r)-sinh(omega*r);
		c1 = a*b*c/(a+b)*i*phi/(i*phi+a*(a+b));
		c2 = sum(exp(-omega*r*radii)./radii)/omega/r;
		c3 = (omega*r+1)*exp(-omega*r);
		null(ii, jj)=1+i*phi-c1+3*c0*c1/(omega*r)^3*(c3-c0*c2);
plot(null); drawnow
	end;
end;

