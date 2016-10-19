function [fxc, fxs, fyc, fys, fzc, fzs, rho, z] = GaussExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, atf, resolution, ring, maxm)
                                                           
% function [fxc, fxs, fyc, fys, fzc, fzs, rho, z] = GaussExc(rhofield,
% zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, atf,
% resolution, ring, maxm) calculates the electric field distribution of a 
% laser beam focused onto a dielectric mirror with characteristics nv dv
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
    maxm = 2;
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

fac = sin(chi).*sqrt(c0);

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
row = ones(size(zv));

if length(focpos)<5
    focpos = [focpos zeros(1,5-length(focpos))];
end
rho0 = sqrt(sum(focpos(2:3).^2));
psi0 = angle(focpos(2)+i*focpos(3));
for j=0:2*maxm
    psi = j/(2*maxm+1)*2*pi; cp = cos(psi); sp = sin(psi);
    %ef = fac.*bessel(0,over*sqrt((fd*n*ss*cp-focpos(4)).^2 + (fd*n*ss*sp-focpos(5)).^2)/max(fd*n*ss)*2.4048).*exp(i*k*ss*rho0*cos(psi-psi0));
    ef = fac.*exp(-om(1)*(fd*n*ss*cp-focpos(4)).^2 - om(2)*(fd*n*ss*sp-focpos(5)).^2 + i*k*ss*rho0*cos(psi-psi0));
    ef = ef.*exp(-i*10*k*fd*n*ss/max(fd*n*ss));
    tmpxt(:,j+1) = (cp^2*tp.*cc+sp^2*ts).*ef;
    tmpxr(:,j+1) = (-cp^2*rp.*cc+sp^2*rs).*ef;
    tmpyt(:,j+1) = cp*sp*(tp.*cc-ts).*ef;
    tmpyr(:,j+1) = cp*sp*(-rp.*cc-rs).*ef;
    tmpzt(:,j+1) = -cp*ss.*tp.*ef;
    tmpzr(:,j+1) = -cp*ss.*rp.*ef;
end
tmpxt = fft(tmpxt,[],2)/(maxm+0.5);
tmpxr = fft(tmpxr,[],2)/(maxm+0.5);
tmpyt = fft(tmpyt,[],2)/(maxm+0.5);
tmpyr = fft(tmpyr,[],2)/(maxm+0.5);
tmpzt = fft(tmpzt,[],2)/(maxm+0.5);
tmpzr = fft(tmpzr,[],2)/(maxm+0.5);

for j=0:maxm
    jj = besselj(j,barg);
    if j==0
        fxc(:,:,1) = jj*((tmpxt(:,1)*row).*ez1 + (tmpxr(:,1)*row).*ez2);
        fyc(:,:,1) = jj*((tmpyt(:,1)*row).*ez1 + (tmpyr(:,1)*row).*ez2);
        fzc(:,:,1) = jj*((tmpzt(:,1)*row).*ez1 + (tmpzr(:,1)*row).*ez2);
    else
        fxc(:,:,j+1) = i^(-j)*jj*(((tmpxt(:,j+1)+tmpxt(:,end-j+1))*row).*ez1 + ((tmpxr(:,j+1)+tmpxr(:,end-j+1))*row).*ez2);
        fyc(:,:,j+1) = i^(-j)*jj*(((tmpyt(:,j+1)+tmpyt(:,end-j+1))*row).*ez1 + ((tmpyr(:,j+1)+tmpyr(:,end-j+1))*row).*ez2);
        fzc(:,:,j+1) = i^(-j)*jj*(((tmpzt(:,j+1)+tmpzt(:,end-j+1))*row).*ez1 + ((tmpzr(:,j+1)+tmpzr(:,end-j+1))*row).*ez2);
        fxs(:,:,j) = i^(-j-1)*jj*(((tmpxt(:,j+1)-tmpxt(:,end-j+1))*row).*ez1 + ((tmpxr(:,j+1)-tmpxr(:,end-j+1))*row).*ez2);
        fys(:,:,j) = i^(-j-1)*jj*(((tmpyt(:,j+1)-tmpyt(:,end-j+1))*row).*ez1 + ((tmpyr(:,j+1)-tmpyr(:,end-j+1))*row).*ez2);
        fzs(:,:,j) = i^(-j-1)*jj*(((tmpzt(:,j+1)-tmpzt(:,end-j+1))*row).*ez1 + ((tmpzr(:,j+1)-tmpzr(:,end-j+1))*row).*ez2);
    end
end

