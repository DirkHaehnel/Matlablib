function [err, c, zz, z] = FlowDiffusionFun(p, t, y, pic, pos)

if length(t)<length(y)
    c = t;
    t = y(:);
    tD = p(1)*t+1;
    zz = [ones(length(t),1) exp(-(t.^2*(p(2)+p(3)))./tD)./tD ...
        exp(-((1-p(2)*t).^2+p(3)*t.^2)./tD)./tD ...
        exp(-((1+p(2)*t).^2+p(3)*t.^2)./tD)./tD ...
        exp(-((1-p(3)*t).^2+p(2)*t.^2)./tD)./tD ...
        exp(-((1+p(3)*t).^2+p(2)*t.^2)./tD)./tD];
    for j=1:size(c,2)
        err(:,j) = zz(:,[1 j+1])*c(:,j);
    end
else
    t = t(:); 
    tD = p(1)*t+1;
    zz = [ones(length(t),1) exp(-(t.^2*(p(2)+p(3)))./tD)./tD ...
        exp(-((1-p(2)*t).^2+p(3)*t.^2)./tD)./tD ...
        exp(-((1+p(2)*t).^2+p(3)*t.^2)./tD)./tD ...
        exp(-((1-p(3)*t).^2+p(2)*t.^2)./tD)./tD ...
        exp(-((1+p(3)*t).^2+p(2)*t.^2)./tD)./tD];

    for j=1:5
        if nargin>4 && ~isempty(pos)
            c(:,j) = lsqnonneg(zz(:,[1 j+1]),y(:,j));
        else
            c(:,j) = zz(:,[1 j+1])\y(:,j);
        end
        z(:,j) = zz(:,[1 j+1])*c(:,j);
    end

    if nargin>3 && ~isempty(pic)
        if pic==1
            subplot(4,1,1:3)
            plot(t, y, 'ob', t, z, 'r'); 
            axis tight
            subplot(4,1,4)
            plot(t, (y-z)./sqrt(abs(z)),t,0*t)
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

