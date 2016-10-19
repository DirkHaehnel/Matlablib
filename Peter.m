% Program for analysis of SMD-TSCPC data

cd c:/joerg/data/peter
% flag = 0; % data-picture (macroscopic versus microscopic time scale)
% flag = 1; % reading data, finding bursts, applying maximum likelihood and plotting figure 1a,b and 2 for the paper
% flag = 2; % calculating the theoretical error probabilities, saving in theory.mat
% flag = 3; % drawing comparison between experiment and theory for rhodamine
% flag = 4; % drawing comparison between experiment and theory for tritc
% flag = 5; % drawing comparison between experiment and theory for mixture
flag = 6; % histogram figure for paper
% flag = 7; % Dick's sorting idea

% names = ['se961100'; 'se961101'; 'se961102']; % rhodamine
% names = ['se961103'; 'se961104'; 'se961105']; % rhodamine
% names = 'se961800'; % tritc
% names = 'se961801'; % tritc
% names = ['se962104'; 'se962105'; 'se962106']; % mixture
taurhod = 19;
tautritc = 12;
dt = 0.1985;

if flag==0
	load se961100;
	[tclen, N] = size(x);
	N = 100;
	tcspc = sum(x');
	data = sum(x);
	[tmp, k1] = max(tcspc(1:(tclen/2)));
	[tmp, k2] = max(tcspc((tclen/2+1):tclen)); k2 = k2 + tclen/2;
	mask = [k1+10:k2-5 k2+10:105];
	t = 1:N;
		
	y = sum(x(mask,t));
	mn = zeros(1,N);
	sig = zeros(1,N);
	lee = 5;
	tmp = zeros(2*lee+1, N-2*lee);
	sig0 = 5;
	for k = 1:lee
		mn(k) = mean(y(1:k+lee));
	end
	for k=1:2*lee+1
		tmp(k,:) = y(k:N-2*lee-1+k);
	end;
	mn(lee+1:N-lee) = mean(tmp);
	for k = N-lee+1:N
		mn(k) = mean(y(k-lee:N));
	end
	for k = 1:lee
	sig(k) = mean((y(1:k+lee)-mn(1:k+lee)).^2);
end
	for k=1:2*lee+1
		tmp(k,:) = (y(k:N-2*lee-1+k) - mn(k:N-2*lee-1+k)).^2;
	end;
	sig(lee+1:N-lee) = mean(tmp);
	for k = N-lee+1:N
		sig(k) = mean((y(k-lee:N)-mn(k-lee:N)).^2);
	end
	
	z = mn + (y-mn).*sig./(sig + sig0);
	tsh = 1.66*mean(z)
	close all
	plot(t, z, [min(t) max(t)], [tsh tsh]);
	tmp = axis;
	axis([min(t) max(t) tmp(3) tmp(4)]);
	xlabel('Time in milliseconds');
	ylabel('Lee filtered data');

	figure
	mygray = ((fix((72:-1:8)/16)/4))'*[1 1 1];
	y = x(:,t);
	mm = max(max(y(mask,:)));
	y(y>mm) = mm*ones(size(y(y>mm)));
	pic = pcolor(t, dt*(1:128), y); set(pic, 'EdgeColor', 'none'); colormap(mygray); axis([0.8 1.001*max(t) dt*0.5 dt*128.2]); h=colorbar; set(h, 'ytick', 0:4); set(h, 'xtick', []); xlabel('Time in milliseconds'); ylabel('Time in nanoseconds'); times
% postscript muss nachbearbeitet werde -> eliminiere nichtganze colorbar-labels!
end

if flag==1 | flag==2
	res = [];
	irf = zeros(128, 1);
	bursts = [];
	burstint = 0;
	irfint = 0;
	for count = 1:length(names(:,1))
		eval(['load ', names(count,:)]);
		[tclen, N] = size(x);
		tcspc = sum(x');
		data = sum(x);
		[tmp, k1] = max(tcspc(1:(tclen/2)));
		[tmp, k2] = max(tcspc((tclen/2+1):tclen)); k2 = k2 + tclen/2;
		mask = [k1+10:k2-5 k2+10:(k1+95)];
		t = 1:N;
		
		y = sum(x(mask,:));
		mn = zeros(1,N);
		sig = zeros(1,N);
		lee = 5;
		tmp = zeros(2*lee+1, N-2*lee);
		sig0 = 5;
		for k = 1:lee
			mn(k) = mean(y(1:k+lee));
		end
		for k=1:2*lee+1
			tmp(k,:) = y(k:N-2*lee-1+k);
		end;
		mn(lee+1:N-lee) = mean(tmp);
		for k = N-lee+1:N
			mn(k) = mean(y(k-lee:N));
		end
		for k = 1:lee
		sig(k) = mean((y(1:k+lee)-mn(1:k+lee)).^2);
		end
		for k=1:2*lee+1
			tmp(k,:) = (y(k:N-2*lee-1+k) - mn(k:N-2*lee-1+k)).^2;
		end;
		sig(lee+1:N-lee) = mean(tmp);
		for k = N-lee+1:N
			sig(k) = mean((y(k-lee:N)-mn(k-lee:N)).^2);
		end
		
		z = mn + (y-mn).*sig./(sig + sig0);
		
		tsh = 1.66*mean(z)
		t1=t([0 z(1:N-1)]<=tsh & z>tsh);
		t2=t([z(2:N) 0]<=tsh & z>tsh);
		
		[tmp, offset] = size(bursts);
		bursts = [bursts zeros(tclen, sum(t2>t1))]; 
		burstint = burstint + sum(t2-t1);
		irfint = irfint + length(t) - sum(t2-t1);

		j = 1;
		for k = 1:length(t1)
			if t2(k)>t1(k)
				bursts(:, offset+j) = sum(x(:, t1(k):t2(k))')';
				j = j+1;
			end
		end
		
		t2 = [1 t2];
		t1 = [t1 N];
		for k = 1:length(t1)
			if t1(k)>t2(k)
				irf = irf + sum(x(:, t2(k):t1(k))')';
			end
		end
		clear x;
	end

	patrhod = convol(irf, 1/taurhod*exp(-(1/taurhod)*((1:k2-k1)-1)));
	patrhod = patrhod'/sum(patrhod(mask));
	pattritc = convol(irf, 1/tautritc*exp(-(1/tautritc)*((1:k2-k1)-1)));
	pattritc = pattritc'/sum(pattritc(mask));
	
	res = [sum(bursts(mask,:)); log(patrhod(mask)./pattritc(mask))'*bursts(mask,:)];
	
	% begin figure
	close all
	t = min(mask):max(mask);
	t1 = mask([mask(2:length(mask)) max(mask)+1]-1~=mask);
	t2 = mask([min(mask)-1 mask(1:length(mask)-1)]+1~=mask);
	semilogy(dt*t, patrhod(t), dt*t, pattritc(t), '--');
	tmp = axis;
	hold
	plot(dt*[t1 t1], [tmp(3) tmp(4)], 'black', dt*[t2 t2], [tmp(3) tmp(4)], 'black');
	hold
%	axis([dt*min(mask) dt*max(mask) tmp(3) tmp(4)]);
	axis([dt*(min(mask)-t1+63) dt*(max(mask)-t1+63) tmp(3) tmp(4)]);
	xlabel('Time in ns'); ylabel('Probability'); times
	text(dt*42, 1.5e-3, 'TRITC');
	text(dt*42, 1e-2, 'Rhodamine 6G');
	figure
%	c = sum(bursts(mask,:)')/patrhod(mask)'/burstint; semilogy(dt*t, irf(t)/irfint, 'x', dt*t, sum(bursts(t,:)')/burstint, 'o', dt*t, c*patrhod(t));
	% Time shifting for congruence reasons
	c = sum(bursts(mask,:)')/pattritc(mask)'/burstint; semilogy(dt*(t-t1+63), irf(t)/irfint, 'x', dt*(t-t1+63), sum(bursts(t,:)')/burstint, 'o', dt*(t-t1+63), c*pattritc (t));
	tmp = axis;
	hold
	plot(dt*[63 63], [tmp(3) tmp(4)], 'black', dt*[t2-t1+63 t2-t1+63], [tmp(3) tmp(4)], 'black');
	hold
	axis([dt*(min(mask)-t1+63) dt*(max(mask)-t1+63) tmp(3) tmp(4)]);
	xlabel('Time in ns'); ylabel('Photoelectron counts per millisecond'); times
end

if flag==2
	step = 0.0025;
	len = 1e3;
	num = 10:10:300;
	s = log(patrhod(mask)./pattritc(mask));
	pa = patrhod(mask); pb = pattritc(mask);
	z = (step*(1:len)).^2;
	clear fa fb
	tic
	for k =1:(len/500);
		tmp = z((k-1)*500+1:k*500);
		fa((k-1)*500+1:k*500) = sum((pa*ones(1,500)).*exp(-i*s*z((k-1)*500+1:k*500)))-1;
		fb((k-1)*500+1:k*500) = sum((pb*ones(1,500)).*exp(-i*s*z((k-1)*500+1:k*500)))-1;
	end
	toc

	M = - 20:0.25:20;
	clear resa resb	
	tic
	for j = 1:length(num)
		for k = 1:length(M)
			resa(j, k) = step*real(sum(sqrt(z).*exp(i*z*M(k)+num(j)*fa)))/2/pi;
			resb(j, k) = step*real(sum(sqrt(z).*exp(i*z*M(k)+num(j)*fb)))/2/pi;
		end
		num(j)
	end
	toc
	save theory M num resa resb patrhod pattritc
end

if flag==3
	load theory
	load rhod
	close all
	subplot(121);
	tmp = [];
	for k=1:length(num)
		for j=1:41
			tmp(k,j) = sum(res(2,res(1,:)>k*10-5 & res(1,:)<k*10+5)>-20+j-1 & res(2,res(1,:)>k*10-5 & res(1,:)<k*10+5)<=-20+j); 
		end
	end
	mygray = flipud(gray);
	h = pcolor(num,min(M):max(M), tmp'); set(h,'EdgeColor', 'none'); 
	hold; plot([min(num) max(num)], [0 0]); hold
	axis([min(num)-2 max(num)+2 1.01*min(M) 1.01*max(M)])
	xlabel('Total number of PE');
	ylabel('ML function value')
	subplot(122);
%	h = pcolor(num,M,resa'); set(h,'EdgeColor', 'none'); colormap(jet);
	h = pcolor(num,M, resa'); set(h,'EdgeColor', 'none'); colormap(mygray);
	hold; plot([min(num) max(num)], [0 0]); hold
	axis([min(num)-2 max(num)+2 1.01*min(M) 1.01*max(M)])
	xlabel('Total number of PE');
	ylabel('ML function value')
	times
	colorbar
	his = hist(res(1,:), num);
	disp(['Exp. Error = ', num2str(sum(res(2,:)>=0)/length(res(2,:)))]);
	disp(['Rhod Error = ', num2str(((sum(resa(:,M>=0)'))*his')/sum(his))]);
end

if flag==4
	load theory
	load tritc
	close all
	mygray = flipud(gray);
	subplot(121);
	tmp = [];
	for k=1:length(num)
		for j=1:41
			tmp(k,j) = sum(res(2,res(1,:)>k*10-5 & res(1,:)<k*10+5)>-20+j-1 & res(2,res(1,:)>k*10-5 & res(1,:)<k*10+5)<=-20+j); 
		end
	end
	h = pcolor(num,min(M):max(M), tmp'); set(h,'EdgeColor', 'none'); 
	hold; plot([min(num) max(num)], [0 0]); hold
	axis([min(num)-2 max(num)+2 1.01*min(M) 1.01*max(M)])
	xlabel('Total number of PE');
	ylabel('ML function value')
	subplot(122);
%	h = pcolor(num,M,resb'); set(h,'EdgeColor', 'none'); colormap(jet);
	h = pcolor(num,M, resb'); set(h,'EdgeColor', 'none'); colormap(mygray);
	hold; plot([min(num) max(num)], [0 0]); hold
	axis([min(num)-2 max(num)+2 1.01*min(M) 1.01*max(M)])
	xlabel('Total number of PE');
	ylabel('ML function value')
	times; colorbar
	his = hist(res(1,:), num);
	disp(['Exp. Error = ', num2str(sum(res(2,:)>=0)/length(res(2,:)))]);
	disp(['Tritc Error = ', num2str(((sum(resa(:,M<=0)'))*his')/sum(his))]);
end

if flag==5
	load theory
	load mixture
	close all
	subplot(121);
%	simple scatter plot
%	plot(res(1,:), res(2,:), '.'); axis([min(num) max(num) min(M) max(M)])
	clear tmp
	for k=1:length(num)
		for j=1:40
			tmp(k,j) = sum(res(2,res(1,:)>k*10-5 & res(1,:)<k*10+5)>-20+j-1 & res(2,res(1,:)>k*10-5 & res(1,:)<k*10+5)<=-20+j); 
		end
	end
	mygray = flipud(gray);
	h = pcolor(num,min(M)+0.5:max(M)-0.5,tmp'); set(h,'EdgeColor', 'none');
	hold; plot([min(num) max(num)], [0 0]); hold
	axis([min(num)-2 max(num)+2 1.01*min(M) 1.01*max(M)])
	xlabel('Total number of PE');
	ylabel('ML function value')
	subplot(122);
	his = hist(res(1,:), num);
	C(1) = sum(res(2,:)>0)/sum(res(1,:)>0);
	C(2) = sum(res(2,:)<0)/sum(res(1,:)>0);
	h = pcolor(num,M, C(1)*resa' + C(2)*resb');
	set(h,'EdgeColor', 'none'); colormap(mygray);
	hold; plot([min(num) max(num)], [0 0]); hold
	axis([min(num)-2 max(num)+2 1.01*min(M) 1.01*max(M)])
	xlabel('Total number of PE');
	ylabel('ML function value')
	times; colorbar
	figure
	[dummy, k] = sort(res(1,:));
	res = fliplr(res(:,k));
	plot(res(1,:),cumsum(res(2,:)>0)./(1:length(res(1,:))))
	times
end

if flag==6
	close all
	load theory
	load rhod
	subplot(311)
	hist(res(2,:), -19.5:19.5);
	hold
	plot(M, 4*sum(resa(round(res(1, res(1,:)<305)/10), : )));
	plot([0 0], [0 100]);
	hold
	axis([-20 20 0 100]);
	text(4,70,'Pure Rhodamine 6G sample');
	disp('rhod error = ')
	sum(res(2,:)<0)/length(res(2,:))
	tmp = sum(resa(round(res(1, res(1,:)<305)/10), : ));
	(sum(tmp(M<0))+0.5*tmp(M==0))/sum(tmp)
	title('        Identified as TRITC          |      Identified as Rhodamine 6G')
	load tritc
	subplot(312)
	hist(res(2,:), -19.5:19.5);
	hold
	plot(M, 4*sum(resb(round(res(1, res(1,:)<305)/10), : )));
	plot([0 0], [0 100]);
	hold
	axis([-20 20 0 100]);
	ylabel('Number of bursts');
	text(4,70,'Pure TRITC sample');
	disp('Tritc error =');
	sum(res(2,:)>0)/length(res(2,:))
	tmp = sum(resb(round(res(1, res(1,:)<305)/10), : ));
	(sum(tmp(M>0))+0.5*tmp(M==0))/sum(tmp)

	load mixture
	subplot(313)
	hist(res(2,:), -19.5:19.5);
	hold
	plot(M, 4*sum(resa(round(res(1,res(2,:)>0)/10), : )) + 4*sum(resb(round(res(1,res(2,:)<0)/10), : )), M, 4*sum(resa(round(res(1,res(2,:)>0)/10), : )), M, 4*sum(resb(round(res(1,res(2,:)<0)/10), : )));
	plot([0 0], [0 100]);
	hold
	axis([-20 20 0 100]);
	xlabel('ML function value'); 
	text(4,70,'Mixture');
	disp('Ratio =');
	sum(res(2,:)>0)
	sum(res(2,:)<0)
	times
	figure
	semilogy(10:10:300, (sum(resa( : ,M<0)')+0.5*resa(: , M==0)'), 10:10:300, (sum(resb( : ,M>0)')+0.5*resb(: , M==0)'));
	text(150, 2e-2, 'TRITC'); text(150, 6e-2, 'Rhodamine 6G');
	xlabel('Total number of photoelectrons per burst'); ylabel('Error rate')
	times
end

if flag==7
	close all
	load theory
	lama = sum(resa(2 ,M<0))+0.5*resa(2 , M==0);
	lamb = sum(resb(2 ,M>0))+0.5*resa(2 , M==0);
	t = (1:10)*20;
a1(1) = 1; b1(1) = 1;
	for k=1:10
		a1(k+1) = a1(k)*(1-lama);
		b1(k+1) = b1(k)*lamb;
	end
	a1(1) = []; b1(1) = [];
	f1 = a1./(a1+b1);
	a2 = 1-(sum(resa(2*(1:10) ,M<0)')+0.5*resa(2*(1:10) , M==0)');
	b2 = sum(resb(2*(1:10) ,M>0)')+0.5*resb(2*(1:10) , M==0)';
	f2 = a2./(a2+b2);
	subplot(211)
	plot(t,f1,'o',t,f2)
	ylabel('purity');
	subplot(212)
	plot(t,a1,'o',t,a2);
	ylabel('amount');
	xlabel('PE equivalent per burst')
end
