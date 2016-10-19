function [err, c, zz, z] = ExpStretchFun(p, t, y, pic)

if length(t)<length(y)
    c = t;
    t = y(:);
    ktau = 1./p(1:end/2); 
    alpha = p(end/2+1:end);
    t = t(:)-min(t); 
    for k=1:length(ktau)
        zz(:,k) = exp(-t.^alpha(k)*ktau(k));
    end
    %zz = zz./(ones(length(t),1)*sum(zz));
    col = ones(length(t),1);%/length(t);
    for j=1:size(c,2)
        err(:,j) = [col zz]*c(:,j);
    end
else
    t = t(:)-min(t); ktau = 1./p(1:end/2); alpha = p(end/2+1:end);
    [m, n] = size(y);
    if m<n y=y'; tmp=n; n=m; m=tmp; end
    for k=1:length(ktau)
        zz(:,k) = exp(-t.^alpha(k)*ktau(k));
    end
    %zz = zz./(ones(m,1)*sum(zz));
    col = ones(length(y),1);%/length(y);

    for j=1:n
        %c(:,j) = lsqnonneg([col zz],y(:,j));
        c(:,j) = [col zz]\y(:,j);
        z(:,j) = [col zz]*c(:,j);
    end

    if nargin>3 && ~isempty(pic)
        if pic==1
            subplot(4,1,1:3)
            plot(t, y, 'ob', t, z, 'r'); 
            axis tight
            subplot(4,1,4)
            plot(t, (y-z)./sqrt(z),t,0*t)
            axis tight
            drawnow
        elseif pic==2
            semilogy(t, y, 'o', t, z); drawnow
        else
            semilogx(t, y, 'o', t, z); drawnow
        end
    end

    err = sum(sum((y-z).^2./abs(z)));
    %err = sum((y-z).^2);
    %ind = y>0;
    %err = sum(y(ind).*log(y(ind)./z(ind))-y(ind)+z(ind));
end

