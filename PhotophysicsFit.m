function [err, z] = PhotophysicsFit(p,int,y,irf,pexc)

kisc = p(1);
kph = p(2);
kcis0 = p(3);
kcis1 = p(4);
ktrans0 = p(5);
ktrans1 = p(6);

col = ones(size(irf,1),1);
row = ones(1,size(y,2));
dia = [1 0; 0 1];
t = cumsum(col);
for j=1:length(int)
    exc = PhotophysicsSat(pexc,int(j)*irf);
    s = [1,0]';
    for k=1:length(exc)-1
        mm = [-exc(k)*(kcis1+kisc)-kcis0-kph exc(k)*ktrans1+ktrans0-kph;...
            exc(k)*kcis1+kcis0 -exc(k)*ktrans1-ktrans0];
        em = expm(mm);
        if kph>0
            s(:,k+1) = em*s(:,k) + (em-dia)*(mm\[kph 0]');
        else
            s(:,k+1) = em*s(:,k);
        end
    end
    bla(:,j) = s(1,:)';
    c = (exc.*s(1,:)')\y(:,j);
    z(:,j) = (exc.*s(1,:)')*c;
    %z(:,j) = cexc*irf.*s(1,:)';
end

if size(y,2)==1
    plot(t,y/max(y),'o',t,z./max(z))
else
    %plot3(t*row,col*(1:size(y,2)),y./(col*max(y)),'.',t*row,ones(size(t))*(1:size(y,2)),z./(col*max(y)))
    subplot(121)
    plot(t,y./(col*max(y)),'o',t,z./(col*max(z)))
    subplot(122)
    plot(t,bla)
end
drawnow

%err = sum(sum((z-y).^2./abs(z)))
err = sum(sum((y./(col*max(y))-z./(col*max(z))).^2));
disp([err p(:)']) 
%err = sum(sum((z-y).^2.))


