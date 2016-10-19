function [err, c, z, zz, back] = AffineFit2D(p, txd, tyd, yd, txm, tym, ym, pic)

% Function AffineFit for affine fitting of a 2D model against data
% (c) Joerg Enderlein, http://www.joerg-enderlein.de (2009)

if length(p)<4
    p(4) = p(3);
end
txm = (txm-p(1))*p(3);
tym = (tym-p(2))*p(4);

zz = interp2(txm,tym,ym,txd,tyd,'spline',0);
back = ones(size(yd));

ind = isfinite(yd);
c = lsqnonneg([back(ind) zz(ind)],yd(ind));
z = reshape([back(:) zz(:)]*c, size(yd,1), size(yd,2));

if nargin>6 && ~isempty(pic)
    %plot(txd,sum(yd),txd,sum(z))
    mim(cat(3,yd,z))
    drawnow
end

err = sum(sum((yd(ind)-z(ind)).^2))

