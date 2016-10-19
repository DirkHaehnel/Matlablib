% Fit of Hong's data
% function param = hongfit(t, z, conc);
load hong
global model
model = 1;
tt=[4:5]; t(tt,:)=[]; z(tt,:)=[]; conc(tt)=[];
t = 60*t;
m = size(t);
ind = sum((t>0)')+1;
close all
p = simplex('hongfun', [1e-5 1e3 0.1], [0 0 0], [], [], [], t, z, conc, ind);
y = zeros(m(1), m(2));
B = p(1); k1 = p(2); ke = p(3)/k1;
for k=1:m(1)
	A = conc(k);
	if model==1
		delta = sqrt((A+B+ke)^2 - 4*A*B)/2;
		x1 = (A+B+ke)/2+delta;
		x2 = (A+B+ke)/2-delta;
		y(k,:) = x1*x2*(exp(2*delta*k1*t(k,:))-1)./(x1*exp(2*delta*k1*t(k,:))- x2);
	end
end
c = reshape(z,1,m(1)*m(2))/reshape(y,1,m(1)*m(2));
plot(t(1,1:ind(1)),c*y(1,1:ind(1)), t',z','o');
hold
for k=2:m(1)
	plot(t(k,1:ind(k)),c*y(k,1:ind(k)));
end
hold
