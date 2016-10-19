function trot = Rad2RotoDiff(rad,tt,vv)

% rad(1) symmetry axis
% rad(2) orthogonal axes

if nargin<2 || isempty(tt)
    tt = 20;
end
if length(tt)<size(rad,1)
    tt = tt*ones(size(rad,1),1);
end
tt(tt<100) = tt(tt<100) + 273.15;
T = 273.15 + (0:10:100);   
BoltzmannConstant = 1.380650324e-23; % NIST
if nargin<3 || isempty(vv)
    visc = [1.7920 1.3080 1.0050 0.8007 0.6560 0.5494 0.4688 0.4061 0.3565 0.3165 0.2838];
    vv = interp1(T,visc,tt,'cubic');
end
if size(rad,2)==1
    trot = 1./(6*tt*BoltzmannConstant/8/pi./vv./rad.^3/1e-30)*1e6; % hydrodyn. radius in Angstrom
else
    d0 = (tt.*BoltzmannConstant./8./pi./vv./rad(:,1)./rad(:,2).^2/1e-30)/1e6; % mean diff coef
    q = rad(:,2)./rad(:,1);
    ind = q<1;
    ft = 0*d0; fra = ft; frb = ft;
    if ~isempty(ind) % prolate ellipsoid, rad(1) long (symmetry) axis, rad(2) short axes
        ft(ind) = sqrt(1-q(ind).^2)./q(ind).^(2/3)./log((1+sqrt(1-q(ind).^2))./q(ind)); 
        fra(ind) = 4*(1-q(ind).^2)./(3*(2-2.*q(ind).^(4/3)./ft(ind))); % symmetry axis
        frb(ind) = 4*(1-q(ind).^4)./q(ind).^2./(3*(2./q(ind).^(2/3).*(2-q(ind).^2)./ft(ind)-2)); % perpendicular axis
    end
    ind = q>1;
    if ~isempty(ind) % oblate elliposid, rad(1) short (symmetry) axis, rad(2) long axes
        ft(ind) = sqrt(q(ind).^2-1)./q(ind).^(2/3)./atan(sqrt(q(ind).^2-1));
        fra(ind) = 4*(1-q(ind).^2)./(3*(2-2*q(ind).^(4/3)./ft(ind)));
        frb(ind) = 4*(1-q(ind).^4)./q(ind).^2./(3*(2./q(ind).^(2/3).*(2-q(ind).^2)./ft(ind)-2));
    end
    fra(q==1) = 1;
    frb(q==1) = 1;
    trot = 1/6*(1./d0*[1 1]).*[fra frb];
end
