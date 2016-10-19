function [err, c, z] = AffineFit2f(p, td, yd, tm, ym);

td = td(:); t = tm(:)*abs(p(:)'); 
if td(1)<min(t(1,:)) 
   t = [td(1)*ones(1,length(p)); t];
   ym = [ones(1,size(ym,2)); ym];
end
if td(end)>max(t(end,:))
   t = [t; td(end)*ones(1,length(p))];
   ym = [ym; zeros(1,size(ym,2))];
end

zz = [];

for j=1:length(p)
    zz = [zz interp1(t(:,j),ym,td,'cubic')];
end

z = 0*yd;
for j=1:size(yd,2)
    c(:,j) = lsqnonneg([ones(size(yd,1),1) zz(:,j)],yd(:,j));
    z(:,j) = [ones(size(yd,1),1) zz(:,j)]*c(:,j);
end
z(:,j) = [ones(size(yd,1),1) zz(:,j)]*c(:,j);

semilogx(td,yd,td,z); drawnow
%plot(td,yd,td,z,t,[ones(length(t),1) ym]*c); drawnow

err = sum(sum((yd-z).^2))

%close; j=1; p=simplex('affinefit',1,0,[],[],[],t,y(:,j),autotime,modres(:,j));
%[err, z] = affinefit(p,t,y(:,j),autotime,modres(:,j));