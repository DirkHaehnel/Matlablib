flag = 2 

if flag==1
	close all
	N = 620;
	e = ones(N,1);
	D = 0.4*spdiags([e -2*e e], -1:1, N, N);
	D(1,1) = D(1,1)- sum(D(1,:));
	D(N,N) = D(N,N) - sum(D(N,:));
	u = zeros(N, 1);
	res = zeros(N, 10);
	u(1:3) = [1 1 1]';
	inert = 150:165;
	for k =1:20000
		tmp = (u-0.2)>0;
		tmp(inert) = zeros(size(inert));
		u = u + 0.001*(tmp-u) + D*u;
		if rem(k,2000)==0
			plot(u); axis([0 N 0 1]); drawnow;
			res(:, k/2000) = u;
		end
	end
	res1 = [ones(1, 10); zeros(inert(1), 10); res(inert, :); zeros(N-inert(length(inert))+1, 10)];
	waterfall(res'); hold; waterfall([1 1:inert(1) inert inert(length(inert)):N], 1:10, res1'); hold; axis([0 N 0 11 0 4]); axis('off'); view([30, 30]);
end

if flag==2
	close all
	t = 0:0.001:1; t = t';
	J = - sqrt(t.^2-2*(t>0.2).*(t-0.2));
	for k=1:7
		tmp = sqrt(t.^2-2*(t>0.2).*(t-0.2)+0.72*exp(0.4*(k-6))-0.12);
		J = [J -tmp tmp];
	end
	J(:, 4:5) = [];
	J0 = sqrt(t.^2-2*(t>0.2).*(t-0.2)+0.72*exp(0.4*(4.5-6))-0.12);
	t0 = 0.2 + (0.72*exp(0.4*(4.5-6))-0.12)/2;
	J0 = t.*(t<t0) + J0.*(t>t0);
	plot(t, real(J), t(t>0.2), -real(J(t>0.2,1)), t, real(J0), '--', [0 0], [-1.2 1.2], [1 1], [-1.2 1.2], [0.2 0.2], [-1.2 1.2], ':'); axis('off');
end
