function [err, c, zz, z] = PulsedExpFun(p, t, y, period, pic)

if length(t)<length(y)
    c = t;
    t = y(:);
    t0 = p(1);
    tau = p(2);
    p = p(3:end);
    p = 1./p(:)';
    t = t(:)-t0; 
    zz = [ones(length(y),1) exp(-t*p)].*(exp(-mod(t,period)/tau)*ones(1,length(p)+1));
    for j=1:size(c,2)
        err(:,j) = zz*c(:,j);
    end
else
    t0 = p(1);
    tau = p(2);
    p = p(3:end);
    t = t(:)-t0; p = 1./p(:)';
    [m, n] = size(y);
    if m<n y=y'; tmp=n; n=m; m=tmp; end
	t = t(isfinite(sum(y,2)));
    y = y(isfinite(sum(y,2)),:);
    [m, n] = size(y);
    zz = [ones(length(y),1) exp(-t*p)].*(exp(-mod(t,period)/tau)*ones(1,length(p)+1));

    for j=1:n
        c(:,j) = lsqnonneg(zz,y(:,j));
        %c(:,j) = [col zz]\y(:,j);
        z(:,j) = zz*c(:,j);
    end

    if nargin>4 && ~isempty(pic)
        %semilogy(t, y, 'o', t, z); drawnow
        plot(t, y, 'ob', t, z, 'r'); drawnow
    end

    %err = sum(sum((y-z).^2.*abs(z)));
    err = sum(sum((y-z).^2));
end

