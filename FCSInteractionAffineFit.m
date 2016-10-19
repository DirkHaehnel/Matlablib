function [err, c, z] = FCSInteractionAffineFit(p, td, yd, tm, ym, pic)

p = p(:); 

td = td(:); t = tm(:)*abs(p(:)'); 
if td(1)<max(t(1,:)) 
   t = [td(1)*ones(1,length(p)); t];
   ym = [ym(1,:); ym];
end
if td(end)>min(t(end,:))
   t = [t; td(end)*ones(1,length(p))];
   ym = [ym; zeros(1,size(ym,2))];
end

for k=1:2
    [tmp, ind] = unique(t(:,k));
    zz(:,k) = interp1(tmp,ym(ind,k),td,'cubic');
end
zz(:,1) = zz(:,1).*zz(:,2);

c = lsqnonneg(zz,yd);
z = zz*c;

if nargin>5 && ~isempty(pic)
    if pic==1
        semilogx(td,yd,td,z);
    else
        plot(td,yd,td,z); 
    end
    drawnow
end

err = sum(sum((yd-z).^2));

