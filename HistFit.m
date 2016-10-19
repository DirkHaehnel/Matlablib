clf;

ind = t(:,6)>1500;
t(ind,:) = [];

ind = t(:,3)>=20;

fe = t(ind, 5) ./ (t(ind, 5) + t(ind, 4));
[h, hx] = hist(fe,0:0.05:1);
bar(hx(hx > 0.1), h(hx>0.1));
xlabel('FRET Efficiency');
ylabel('# Bursts');

s = ['# donor only bursts (E <= 0.1) = ' int2str(sum(h(hx<=0.1)))];
text(0.1,0.95,s,'VerticalAlignment','Top','units','normal');

figure;
p0 = [0.2 0.7 0.1 0.1];
p = Simplex('Gauss',p0,[],[],[],[],hx(hx > 0.35), h(hx>0.35),[],[],1)

%p0 = [0.1 0.45 0.1 0.1];
%p = Simplex('Gauss',p0,[],[],[],[],hx(hx > 0.1 & hx < 0.7), h(hx>0.1 & hx < 0.7),[],[],1)