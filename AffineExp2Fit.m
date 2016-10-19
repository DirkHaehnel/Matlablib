function [err, c, z] = AffineExp2Fit(p, td, yd, tm, ym);

td = td(:); tm = p(1)*tm(:); p=p(:)';
if td(1)<tm(1) 
   tm = [td(1); tm];
   ym = [ones(1,size(ym,2)); ym];
end
if td(end)>tm(end)
   tm = [tm; td(end)];
   ym = [ym; zeros(1,size(ym,2))];
end

z = interp1(tm,ym,td,'cubic');
for j=1:size(yd,2)
    c(:,j) = lsqnonneg([ones(size(yd,1),1) exp(-td*p(2)).*z(:,j) z(:,j) exp(-td*p(3))],yd(:,j));
    z(:,j) = [ones(size(yd,1),1) exp(-td*p(2)).*z(:,j) z(:,j) exp(-td*p(3))]*c(:,j);
end
semilogx(td,yd,td,z); drawnow

err = sum(sum((yd-z).^2))

%close; j=1; p=simplex('affinefit',1,0,[],[],[],t,y(:,j),autotime,modres(:,j));
%[err, z] = affinefit(p,t,y(:,j),autotime,modres(:,j));