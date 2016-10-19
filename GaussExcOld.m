function [fxc, fxs, fyc, fys, fzc, fzs, rho, z] = GaussExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, atf, resolution, ring, maxm)
                                                           
% function [exc, rho, z, fx0, fx2, fz] = NewExc(rhofield, zfield, NA, n0, n, n1, d0, d, d1, lamex, over, focpos, atf)
% calculates the electric field distribution of a laser beam focused onto a dielectric mirror with characteristics nv dv
%
% input arguments:
% rhofield - two-element vector containing the minimum and maximum value of the rho coordinate
% zfield   - two-element vector containing the minimum and maximum value of the z coordinate
% NA       - numerical aperture
% fd       - focal distance
% n0       - refractive index vector of layers below sample
% n        - refractive index of sample medium
% n1       - refractive index vector of layers above sample
% d0       - thickness vector of layers below sample
% d        - thickness of sample medium
% d1       - thickness vector of layers above sample
% lamex    - excitation wavelength
% over     - laser characteristics: [beam waist position, beam waist radius, astigmatism value]
% focpos   - focus position [z, x displacement, y displacment, x inclination, y inclination]
% atf      - cover slide correction
%
% output arguments:
% exc      - excitation intensity distribution
% rho, z   - rho and z coordinate grid
% f**      - electric field component distributions
%
% field is given by
% ex = fxc(:,:,1) + fxc(:,:,2)*cos(phi) + ... + fxs(:,:,1)*sin(phi) + ...
% etc.

maxnum = 1e3;
if nargin < 17 | isempty(maxm)
    maxm = 3;
end

if nargin<15 | isempty(resolution)
    resolution = [20 20];
end
if length(resolution)==1
    resolution = resolution*[1 1];
end

if isempty(d) | d<max(zfield(:))
    d = max(zfield(:));
end

n0 = n0(:).'; n1 = n1(:).'; d0 = d0(:).'; d1 = d1(:).';

k0 = 2*pi/lamex*n0; k = 2*pi/lamex*n; k1 = 2*pi/lamex*n1; 
d0 = 2*pi*d0/lamex; d1 = 2*pi*d1/lamex;

if length(rhofield)==1
    dz = lamex/resolution(2); 
    rhov = rhofield; 
    if length(zfield)>1
        zv = zfield(1) + (0.5:(zfield(2)-zfield(1))/dz)*dz;
    else
        zv = zfield;
    end
    if isempty(rhov)
        rhov = rhofield;
    end
    if isempty(zv)
        zv = zfield;
    end
    rho = rhov;
    z = zv;
elseif length(rhofield)==2
    drho = lamex/resolution(1);
    dz = lamex/resolution(2);
    rhov = rhofield(1) + (0.5:(rhofield(2)-rhofield(1))/drho)*drho; 
    if length(zfield)>1
        zv = zfield(1) + (0.5:(zfield(2)-zfield(1))/dz)*dz;
    else
        zv = zfield;
    end
    if isempty(rhov)
        rhov = rhofield;
    end
    if isempty(zv)
        zv = zfield;
    end
    [rho, z] = ndgrid(rhov, zv);
else
    drho = diff(rhofield(1:2,1)); dz = diff(zfield(1,1:2));
    rhov = rhofield(:,1)'; zv = zfield(1,:);
    rho = rhofield; z = zfield;
end

chimax = asin(NA/n0(1)); 
if nargin<15 | isempty(ring) 
    chimin = 0;
else
    chimin = asin(ring/n0(1));    
end

dchi = (chimax-chimin)/maxnum;
chi = (chimin:dchi:chimax)';

c0 = cos(chi);
s0 = sin(chi);
ss = n0(1)/n*s0;
cc = sqrt(1-ss.^2);
rad = NA*fd*s0/sin(chimax); % relation between mode angle and entry ray lateral distance

if length(over)==2 % over = [z0, ww]
    if over(2)==0 % wide field illumination - not an elegant solution but it works
        fac = 0*chi;    
        fac(1) = 1;
        zeta = 0;
    else
        if over(1)==inf
            zeta = 0;
            ww = over(2);
        else
            zeta = over(1)*lamex/pi/over(2)^2;
            ww = over(2)*sqrt(1+zeta^2);
        end
        fac = sin(chi).*sqrt(c0).*exp(-rad.^2/ww^2*(1-i*zeta));
    end
    om = [];
    zeta = [1 1]*zeta;
    w1 = ww; w2 = ww;
elseif length(over)==3 % over = [ast, w0, 1]
    zeta = 2*over(1)*pi*over(2)^2/(NA*fd)^2;
    w0 = over(2)/sqrt(1+zeta^2);
    zeta = [-zeta zeta];
    om = (1-i*zeta)./(1+zeta.^2)./w0^2/2;
    if length(focpos)==1
        fac = sin(chi).*sqrt(c0).*exp(-rad.^2*sum(om));
    else
        fac = sin(chi).*sqrt(c0).*exp(-rad.^2*sum(om)-2*sum(focpos(2:3).^2.*om));
    end
    om = -i*diff(om)*rad.^2;
    w1 = w0; w2 = w0;
else % over = [z1, z2, w1, w2]
    zeta = over(1:2)*lamex/pi./over(3:4).^2;
    om = (1-i*zeta)./(1+zeta.^2)./over(3:4).^2/2;
    if length(focpos)==1
        fac = sin(chi).*sqrt(c0).*exp(-rad.^2*sum(om));
    else
        fac = sin(chi).*sqrt(c0).*exp(-rad.^2*sum(om)-2*sum(focpos(2:3).^2.*om));
    end
    om = -i*diff(om)*rad.^2;
    w1 = over(3); w2 = over(4);
end

[rp0,rs0,tp0,ts0] = Fresnel(n0(1)*c0,[n0 n],d0);
[rp0,rs0] = Fresnel(n*cc,[n n0(end:-1:1)],d0(end:-1:1));
[rp1,rs1] = Fresnel(n*cc,[n n1],d1);
if nargin>13 & ~isempty(atf)
    if length(atf)==1 % accounting for reflection losses when using water immersion; atf = ref index of cover slide
        [tmp, tms, tmp, tms] = Fresnel(n0(1)*c0,n0(1),atf);
        [mp, ms, mp, ms] = Fresnel(sqrt(atf^2-n0(1)^2*s0.^2),atf,n0(1));
    else % + aberration for water immersion; atf(2) = thickness mismatch
        [tmp, tms, tmp, tms] = Fresnel(n0(1)*c0,[n0(1) n0(1) atf(1)], -2*pi*atf(2)/lamex);        
        [mp, ms, mp, ms] = Fresnel(sqrt(atf(1)^2-n0(1)^2*s0.^2),[atf(1) atf(1) n], 2*pi*atf(2)/lamex);        
    end
    tp0 = tmp.*mp.*tp0;
    ts0 = tms.*ms.*ts0;
end
tp = tp0./(1-rp1.*rp0.*exp(2*i*k*cc.'*d));
ts = ts0./(1-rs1.*rs0.*exp(2*i*k*cc.'*d));
tp = tp.';
ts = ts.';
rp = tp.*rp1.';
rs = ts.*rs1.';

phase1 = k*cc*zv - k0(1)*c0*focpos(1)*ones(size(zv));
phase2 = k*cc*(2*d-zv) - k0(1)*c0*focpos(1)*ones(size(zv));
ez1 = exp(i*phase1);
ez2 = exp(i*phase2);

barg = k*rhov'*ss';

if length(over)==2 & length(focpos)==1
    j0 = besselj(0,barg); j1 = besselj(1,barg); j2 = besselj(2,barg);
    fxc = j0*(((fac.*(tp.*cc+ts))*ones(size(zv))).*ez1 + ((fac.*(-rp.*cc+rs))*ones(size(zv))).*ez2); % cos(0*phi)-component
    fxc(:,:,2) = 0;
    fxc(:,:,3) = j2*(-((fac.*(tp.*cc-ts))*ones(size(zv))).*ez1 + ((fac.*(rp.*cc+rs))*ones(size(zv))).*ez2); % cos(2*phi)-component
    fys(:,:,2) = fxc(:,:,3);
    fxs = 0*fys;
    fyc = 0*fxc;
    fzc = 0*fxc;
    fzc(:,:,2) = 2*i*j1*(((fac.*tp.*ss)*ones(size(zv))).*ez1 + ((fac.*rp.*ss)*ones(size(zv))).*ez2); % cos(phi)-component
    fzs = 0*fys;
else
    % Inclination
    if length(focpos)<4
        ac = ones(length(rad),1);
        iac = 0;
        as = [];
        ias = [];
    else
        if focpos(5)==0
            kappa = 2*rad*focpos(4)*(zeta(1)+i)/(1+zeta(1)^2)/w1^2;
            ac(:,1) = besselj(0,kappa);
            for j=1:maxm
                ac(:,j+1) = 2/i^j*besselj(j,kappa);
            end
            as = [];
            iac = 0:maxm;
            ias = [];
        elseif focpos(4)==0
            kappa = 2*rad*focpos(5)*(zeta(2)+i)/(1+zeta(2)^2)/w2^2;
            ac(:,1) = besselj(0,kappa);
            for j=1:maxm
                if mod(j,2)==0
                    ac(:,j/2+1) = 2*besselj(j,kappa);
                else
                    as(:,(j+1)/2) = -2*i*besselj(j,kappa);
                end
            end
            ias = 1:2:maxm;
            iac = 0:2:maxm;
        else
            kappa = 2*rad*focpos(4)*(zeta(1)+i)/(1+zeta(1)^2)/w1^2;
            ac(:,1) = besselj(0,kappa);
            for j=1:maxm
                ac(:,j+1) = 2/i^j*besselj(j,kappa);
            end
            as = [];
            iac = 0:maxm;
            ias = [];
            kappa = 2*rad*focpos(5)*(zeta(2)+i)/(1+zeta(2)^2)/w2^2;
            bc(:,1) = besselj(0,kappa);
            for j=1:maxm
                if mod(j,2)==0
                    bc(:,j/2) = 2*besselj(j,kappa);
                else
                    bs(:,(j+1)/2) = -2*i*besselj(j,kappa);
                end
            end
            ibs = 1:2:maxnum;
            ibc = 0:2:maxm;
            [iac, ias, ac, as] = TrigTimes(iac, ias, ibc, ibs, ac, as, bc, bs, maxm, maxm);
        end
    end
    
    % Astigmatism
    if ~isempty(om) | om==0
        clear bc bs ibc ibs
        bc(:,1) = besselj(0,om);
        for j=1:maxm
            bc(:,j+1) = 2/i^j*besselj(j,om);
        end
        ibc = 0:2:2*maxm;
        bs = []; ibs = [];
        [iac, ias, ac, as] = TrigTimes(iac, ias, ibc, ibs, ac, as, bc, bs, maxm, maxm);
    end
    
    % Displacement
    if length(focpos)>1
        clear bc bs ibc ibs
        om = 2*pi/lamex*rad*sqrt(sum(focpos(2:3).^2))/fd;
        phi = angle(focpos(2)+i*focpos(3));
        bc(:,1) = 0.5*besselj(0,om);
        for j=1:maxm
            tmp = besselj(j,om);
            bc(:,j+1) = i^j*cos(j*phi)*tmp;
            bs(:,j) = i^j*sin(j*phi)*tmp;
        end
        ibc = 0:maxm;
        ibs = 1:maxm;
        [iac, ias, ac, as] = TrigTimes(iac, ias, ibc, ibs, ac, as, bc, bs, maxm, maxm);
    end
        
    % integration over psi
    fac = fac/2; % for achieving uniform normalization across all calculation modes
    [ixc1, ixs1, bc, bs] = TrigTimes(iac, ias, [0 2], [], ac, as, [fac.*(tp.*cc+ts) fac.*(tp.*cc-ts)], [], maxm, maxm);
    xc1 = permute(repmat(bc,[1 1 length(rhov)]),[3 1 2]);
    xs1 = permute(repmat(bs,[1 1 length(rhov)]),[3 1 2]);
    [ixc2, ixs2, bc, bs] = TrigTimes(iac, ias, [0 2], [], ac, as, [fac.*(-rp.*cc+rs) -fac.*(rp.*cc+rs)], [], maxm, maxm);
    xc2 = permute(repmat(bc,[1 1 length(rhov)]),[3 1 2]);
    xs2 = permute(repmat(bs,[1 1 length(rhov)]),[3 1 2]);
    [iyc1, iys1, bc, bs] = TrigTimes(iac, ias, [], 2, ac, as, [], fac.*(tp.*cc-ts), maxm, maxm);
    yc1 = permute(repmat(bc,[1 1 length(rhov)]),[3 1 2]);
    ys1 = permute(repmat(bs,[1 1 length(rhov)]),[3 1 2]);
    [iyc2, iys2, bc, bs] = TrigTimes(iac, ias, [], 2, ac, as, [], -fac.*(rp.*cc+rs), maxm, maxm);
    yc2 = permute(repmat(bc,[1 1 length(rhov)]),[3 1 2]);
    ys2 = permute(repmat(bs,[1 1 length(rhov)]),[3 1 2]);
    [izc1, izs1, bc, bs] = TrigTimes(iac, ias, 1, [], ac, as, 2*fac.*tp.*ss, [], maxm, maxm);
    zc1 = permute(repmat(bc,[1 1 length(rhov)]),[3 1 2]);
    zs1 = permute(repmat(bs,[1 1 length(rhov)]),[3 1 2]);
    [izc2, izs2, bc, bs] = TrigTimes(iac, ias, 1, [], ac, as, 2*fac.*rp.*ss, [], maxm, maxm);
    zc2 = permute(repmat(bc,[1 1 length(rhov)]),[3 1 2]);
    zs2 = permute(repmat(bs,[1 1 length(rhov)]),[3 1 2]);

    % convolution with radial dependency
    bb = besselj(0,barg);
    fxc = zeros(size(rho,1),size(z,2),maxm+1); 
    fxs = zeros(size(rho,1),size(z,2),maxm); 
    fyc = fxc; fys = fxs; fzc = fxc; fzs = fxs;
    if ismember(0,ixc1) fxc(:,:,1) = 2*(bb.*xc1(:,:,ixc1==0))*ez1; end
    if ismember(0,ixc2) fxc(:,:,1) = fxc(:,:,1) + 2*(bb.*xc2(:,:,ixc2==0))*ez2; end
    if ismember(0,iyc1) fyc(:,:,1) = 2*(bb.*yc1(:,:,iyc1==0))*ez1; end
    if ismember(0,iyc2) fyc(:,:,1) = fyc(:,:,1) + 2*(bb.*yc2(:,:,iyc2==0))*ez2; end
    if ismember(0,izc1) fzc(:,:,1) = 2*(bb.*zc1(:,:,izc1==0))*ez1; end
    if ismember(0,izc2) fzc(:,:,1) = fzc(:,:,1) + 2*(bb.*zc2(:,:,izc2==0))*ez2; end
    for j=1:maxm
        bb = 2/i^j*besselj(j,barg);
        
        if ismember(j,ixc1) fxc(:,:,j+1) = (bb.*xc1(:,:,ixc1==j))*ez1; end
        if ismember(j,ixc2) fxc(:,:,j+1) = fxc(:,:,j+1) + (bb.*xc2(:,:,ixc2==j))*ez2; end
        if ismember(j,iyc1) fyc(:,:,j+1) = (bb.*yc1(:,:,iyc1==j))*ez1; end
        if ismember(j,iyc2) fyc(:,:,j+1) = fyc(:,:,j+1) + (bb.*yc2(:,:,iyc2==j))*ez2; end
        if ismember(j,izc1) fzc(:,:,j+1) = (bb.*zc1(:,:,izc1==j))*ez1; end
        if ismember(j,izc2) fzc(:,:,j+1) = fzc(:,:,j+1) + (bb.*zc2(:,:,izc2==j))*ez2; end
        
        if ismember(j,ixs1) fxs(:,:,j) = (bb.*xs1(:,:,ixs1==j))*ez1; end
        if ismember(j,ixs2) fxs(:,:,j) = fxs(:,:,j) + (bb.*xs2(:,:,ixs2==j))*ez2; end
        if ismember(j,iys1) fys(:,:,j) = (bb.*ys1(:,:,iys1==j))*ez1; end
        if ismember(j,iys2) fys(:,:,j) = fys(:,:,j) + (bb.*ys2(:,:,iys2==j))*ez2; end
        if ismember(j,izs1) fzs(:,:,j) = (bb.*zs1(:,:,izs1==j))*ez1; end
        if ismember(j,izs2) fzs(:,:,j) = fzs(:,:,j) + (bb.*zs1(:,:,izs1==j))*ez2; end

    end
             
end
