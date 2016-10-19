function [err, c, zz] = GaussBen(p, z, bck)

[x,y] = meshgrid(-(size(z,2)-1)/2:(size(z,2)-1)/2,-(size(z,1)-1)/2:(size(z,1)-1)/2);

mx = p(1); my = p(2); p(1:2)=[];
wx = p(1:2:end); wy = p(2:2:end);

M = ones(numel(z),1);
if nargin>2 && ~isempty(bck)
    if bck == 1
        M = [M x(:) y(:)];
    end
    if bck == 2
        M = [M x(:) y(:) x(:).^2 y(:).^2 x(:).*y(:)];
    end
end
for j=1:length(wx)
    tmp = exp(-(x-mx).^2/2/wx(j)^2 - (y-my).^2/2/wy(j)^2);
    M = [M tmp(:)];
end
%c = lsqnonneg(M(isfinite(z),:),z(isfinite(z)));
c = M(isfinite(z),:)\z(isfinite(z));
zz = reshape(M*c,size(z,1),size(z,2));
surf(x,y,(z-zz)./abs(zz)); drawnow;

err = sum(sum((z(isfinite(z))-zz(isfinite(z))).^2./abs(zz(isfinite(z)))))

