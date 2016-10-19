function [int, x, y, feldx, feldy, feldz] = DefExcImage(al,be,nn,pixel,rho,fxc,fxs,fyc,fys,fzc,fzs,psi)

if length(nn)==1
    nn = [nn nn];
end

[x,y] = meshgrid(-nn(2):nn(2),-nn(1):nn(1));
p = angle(x+1i*y);
if nargin>11 && ~isempty(psi)
    p = p - psi;
end
r = sqrt(x.^2+y.^2);
rho = rho(:,1)/pixel;

feldx = zeros(size(x,1),size(x,2),size(fxc,2));
feldy = feldx; feldz = feldx;
for k=1:size(fxc,2)
    feldx(:,:,k) = interp1(rho,squeeze(fxc(:,k,1)),r,'spline',0);
    feldy(:,:,k) = interp1(rho,squeeze(fyc(:,k,1)),r,'spline',0);
    feldz(:,:,k) = interp1(rho,squeeze(fzc(:,k,1)),r,'spline',0);
    for j=1:size(fxs,3)
        feldx(:,:,k) = feldx(:,:,k) + interp1(rho,squeeze(fxc(:,k,j+1)),r,'spline',0).*cos(j*p);
        feldx(:,:,k) = feldx(:,:,k) + interp1(rho,squeeze(fxs(:,k,j)),r,'spline',0).*sin(j*p);
        feldy(:,:,k) = feldy(:,:,k) + interp1(rho,squeeze(fyc(:,k,j+1)),r,'spline',0).*cos(j*p);
        feldy(:,:,k) = feldy(:,:,k) + interp1(rho,squeeze(fys(:,k,j)),r,'spline',0).*sin(j*p);
        feldz(:,:,k) = feldz(:,:,k) + interp1(rho,squeeze(fzc(:,k,j+1)),r,'spline',0).*cos(j*p);
        feldz(:,:,k) = feldz(:,:,k) + interp1(rho,squeeze(fzs(:,k,j)),r,'spline',0).*sin(j*p);
    end
end
feldx(isnan(feldx)) = 0;
feldy(isnan(feldy)) = 0;
feldz(isnan(feldz)) = 0;

int = abs(feldx*sin(al)*cos(be) + feldy*sin(al)*sin(be) + feldz*cos(al)).^2;
