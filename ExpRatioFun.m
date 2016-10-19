function [err, c, zz, z] = ExpRatioFun(p, t, y, pic, pos)

if length(t)<length(y)
    c = t;
    t = y(:);
    p = p(:)';
    k = 1./p(1:end/2);
    al = p(end/2+1:end);
    t = t(:)-min(abs(t(:))); 
    zz = [ones(m,1) 1./(1 + (ones(length(t),1)*al).*exp(-abs(t)*k))];
    for j=1:size(c,2)
        err(:,j) = zz*c(:,j);
    end
else
    t = t(:)-min(abs(t(:))); 
    p = p(:)';
    k = 1./p(1:end/2);
    al = p(end/2+1:end);
    [m, n] = size(y);
    if m<n y=y'; tmp=n; n=m; m=tmp; end
	t = t(isfinite(sum(y,2)));
    y = y(isfinite(sum(y,2)),:);
    [m, n] = size(y);
    zz = [ones(m,1) 1./(1 + (ones(length(t),1)*al).*exp(-abs(t)*k))];

    for j=1:n
        if nargin>4 && ~isempty(pos)
            c(:,j) = lsqnonneg(zz,y(:,j));
        else
            c(:,j) = zz\y(:,j);
        end
        z(:,j) = zz*c(:,j);
    end

    if nargin>3 && ~isempty(pic)
        if pic==1
            plot(t, y, 'ob', t, z, 'r'); drawnow
        elseif pic==2
            semilogy(t, y, 'o', t, z); drawnow
        else
            semilogx(t, y, 'o', t, z); drawnow
        end
    end

    err = sum(sum((y-z).^2./abs(z)));
end

