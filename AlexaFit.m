function [err, c, z] = AlexaFit(p, exc, td, yd, tm, ym);

td = td(:); tm = p(1)*tm(:); exc = exc(:);

if td(1)<tm(1) 
   tm = [td(1); tm];
   ym = [ones(1,size(ym,2)); ym];       
end
if td(end)>tm(end)
   tm = [tm; td(end)];
   ym = [ym; zeros(1,size(ym,2))];
end

if exc(4)>0
    exc(6) = ((exc(4)+2*exc(5)+exc(6))./exc(4)-((exc(4)+exc(5))./exc(4)).^2)/2;
    exc(5) = (exc(4)+exc(5))./exc(4);
    exc(3:4) = exc(3:4)./sum(exc(3:4));
else
    exc(3) = 1;
    exc(4:6) = 0;
end

z = interp1(tm,ym,td,'cubic');
for j=1:size(yd,2)
    %vec = [ones(size(yd,1),1) z(:,j) exp(-p(1+j)*td)];
    vec = [ones(size(yd,1),1) z(:,j) exp(-p(1+j)*td).*z(:,j)];
    %vec = [ones(size(yd,1),1) z(:,j) exp(-p(1+j)*exc(5)*td-p(1+j)^2*exc(6)*td.^2)];
    %vec = [ones(size(yd,1),1) (exc(3)+exc(4)*exp(-p(1+j)*exc(5)*td-p(1+j)^2*exc(6)*td.^2)).*z(:,j)];
    c(:,j) = lsqnonneg(vec,yd(:,j));
    z(:,j) = vec*c(:,j);
end
semilogx(td,yd,td,z); drawnow

err = sum(sum(((yd-z).^2)./abs(z)))

