function [err, c, c0, zz, z0, z] = Antibunch(p, t, x, period, bck, bld)

t = t(:); 
if size(x,2)>1
    t = [-t(end:-1:2); t];
    x = [x(end:-1:2,2:2:end); x(:,1:2:end)];
elseif size(x,2)==4
    t = [-t(end:-1:2); t];
    x = x(:,1:2)-x(:,3:4);
    x = [x(end:-1:2,2); x(:,1)];
end
p = p(:)';
if nargin>4 && ~isempty(bck)
    fac = p(end-2); pos = p(end-1); sig = p(end); p(end-2:end) = []; 
end
p(2:end) = 1./p(2:end);

t = t-p(1);
col = ones(length(x),1)/length(x);
if nargin<4 
    zz = -exp(-abs(t)*p(2:end));
    zz = zz./(ones(size(zz,1),1)*sum(zz));
else
    zz = [exp(-abs(mod(t,period))*p(2:end)) + exp((abs(mod(t,period))-period)*p(2:end))];
    znorm = sum(zz);
    zz = [zz -exp(-abs(t)*p(2:end))];
    zz = zz./(ones(size(zz,1),1)*[znorm znorm]);
    if size(x,2)==2
        z0 = zz(:,1:end/2);
    end
end

if nargin>4 && ~isempty(bck)
    zz = [zz abs(t).^fac.*exp(-abs(t-pos)/sig) abs(t).^fac.*exp(-abs(t+pos)/sig)]; 
end

if size(x,2)==2
    %     c = ([col zz]\x(:,1));
    %     c0 = ([col z0]\x(:,2));
    c = ([[col;col] [zz; zeros(size(zz))] [zeros(size(z0)); z0]]\[x(:,1); x(:,2)]);
    c0 = c([1 end-size(z0,2)+1:end]);
    c = c([1:end-size(z0,2)]);
    z = [col zz]*c;
    z = [z [col z0]*c0];
else
    c = ([col zz]\x);
    z = [col zz]*c;
    c0 = 0;
end
if nargin>4 && ~isempty(bck)
    c(end-(2*bck-1):end) = [];
end

if nargin>5 && ~isempty(bld)
    plot(t, x, 'o', t, z, 'markersize', 2.5); drawnow
end

err = sum(sum((x-z).^2./abs(z)));
%err = sum(sum((x-z).^2))



