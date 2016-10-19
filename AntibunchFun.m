function [err, c, c0, zz, z0, z] = AntibunchFun(p, t, x, ctau, tau, period, bck, bld);

ctau = ctau(:);
tau = tau(:)';
t = t(:); 
if size(x,2)>1
    t = [-t(end:-1:2); t];
    x = [x(end:-1:2,2:2:end); x(:,1:2:end)];
end
p = p(:)';
if nargin>6 & ~isempty(bck)
    fac = p(end-2); pos = p(end-1); sig = p(end); p(end-2:end) = []; 
end
p(2:end) = 1./p(2:end);

t = t-p(1);
col = ones(length(x),1)/length(x);
if nargin<6
    zz = -exp(-abs(t)*(1./tau))*ctau;
    zz = zz./(ones(size(zz,1),1)*sum(zz));
else
    zz = (exp(-abs(mod(t,period))*(1./tau)) + exp((abs(mod(t,period))-period)*(1./tau)))*ctau;
    znorm = sum(zz);
    zz = [zz -exp(-abs(t)*(1./tau))*ctau];
    zz = zz./(ones(size(zz,1),1)*[znorm znorm]);
    if size(x,2)==2
        z0 = zz(:,1:end/2);
    end
end

if nargin>6 & ~isempty(bck)
    zz = [zz abs(t).^fac.*exp(-abs(t-pos)/sig) abs(t).^fac.*exp(-abs(t+pos)/sig)]; 
end

if size(x,2)==2
    c = ([col zz]\x(:,1));
    z = [col zz]*c;
    c0 = ([col z0]\x(:,2));
    z = [z [col z0]*c0];
else
    c = ([col zz]\x);
    z = [col zz]*c;
    c0 = 0;
end
if nargin>6 & ~isempty(bck)
    c(end-(2*bck-1):end) = [];
end

if nargin>7 & ~isempty(bld)
    plot(t, x, 'o', t, z, 'markersize', 2.5); drawnow
end

err = sum(sum((x-z).^2./abs(z)));
%err = sum(sum((x-z).^2))



