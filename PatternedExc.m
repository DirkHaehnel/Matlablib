function [fx, fy, fz, exc, tx, ty, tz] = PatternedExc(x, y, pattern, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, incoherent, circular)

dd = diff(x(1,1:2));
if length(zfield)==3
    resolution = [lamex/dd lamex/zfield(3)];
    zfield = [zfield(1)-lamex/resolution(2)/2 zfield(2)];
elseif length(zfield)==2
    resolution = [lamex/dd lamex/0.5];
    zfield = [zfield(1)-lamex/resolution(2)/2 zfield(2)];
else
    resolution = lamex/dd;
end
rhofield = [-lamex/resolution(1)/2 max(sqrt(x(:).^2+y(:).^2))];

over = inf;
focpos = 0;
atf = []; 
exc = GaussExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, atf, resolution);
if nargin>14 && ~isempty(circular) 
    tmp = RotateEMField(exc,pi/2);
    exc.fxc = exc.fxc + i*tmp.fxc;
    exc.fxs = exc.fxs + i*tmp.fxs;
    exc.fyc = exc.fyc + i*tmp.fyc;
    exc.fys = exc.fys + i*tmp.fys;
    exc.fzc = exc.fzc + i*tmp.fzc;
    exc.fzs = exc.fzs + i*tmp.fzs;
end
[tx, ty, tz] = GaussExc2Grid(x,y,exc);
if nargin<14 || isempty(incoherent) || incoherent==0
    fx = fft2(ifftshift(ifftshift(tx,1),2));
    fy = fft2(ifftshift(ifftshift(ty,1),2));
    fz = fft2(ifftshift(ifftshift(tz,1),2));
else
    fx = fft2(ifftshift(ifftshift(abs(tx).^2,1),2));
    fy = fft2(ifftshift(ifftshift(abs(ty).^2,1),2));
    fz = fft2(ifftshift(ifftshift(abs(tz).^2,1),2));
end

ef = repmat(fft2(pattern),[1 1 size(exc.z,2)]);
fx = ifft2(ef.*fx);
fy = ifft2(ef.*fy);
fz = ifft2(ef.*fz);
if nargin>13 && ~(isempty(incoherent) || incoherent==0)
    fx = real(fx);
    fy = real(fy);
    fz = real(fz);
end

