% program SPM-EM
close all

rhofield = [0 2];
height = 3;
zfield = [-height height];
NA = 1.2;
theta_min = 0*asin(0.9*1.2/1.33); % Leuchs
fd = 3e3;
n0 = 1.2;
n = 1.33;
n1 = 1.33;
d0 = [];
d = height;
d1 = [];
lamex = 0.63;
over = 1e6;
focpos = 0;
atf = [];
resolution = 30;

eflag = true;
% eflag = false;

if eflag % electric field
    if 0
        ring = ['cos(psi).*(cos(psi*3)<0.9).*(rad>sin(' num2str(theta_min) '))']; maxm = 50;
        [fxc, fxs, fyc, fys, fzc, fzs, rho, z] = GaussExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, atf, resolution, ring, maxm);
        ring = ['cos(psi).*(cos((pi/2+psi)*3)<0.9).*(rad>sin(' num2str(theta_min) '))'];
        [fxc2, fxs2, fyc2, fys2, fzc2, fzs2] = GaussExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, atf, resolution, ring, maxm);
        [fxc2, fxs2, fyc2, fys2, fzc2, fzs2] = RotateEMField(fxc2, fxs2, fyc2, fys2, fzc2, fzs2, pi/2);
    else
        ring = ['cos(psi).*(rad>sin(' num2str(theta_min) '))']; maxm = 3;
        [fxc, fxs, fyc, fys, fzc, fzs, rho, z] = GaussExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, atf, resolution, ring, maxm);
        [fxc2, fxs2, fyc2, fys2, fzc2, fzs2] = RotateEMField(fxc, fxs, fyc, fys, fzc, fzs, pi/2);
    end
else % magnetic field
    if 0
        ring = ['-sin(psi).*(cos((pi/2+psi)*3)<0.9).*(rad>sin(' num2str(theta_min) '))']; maxm = 50;
        [fxc, fxs, fyc, fys, fzc, fzs, rho, z] = GaussExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, atf, resolution, ring, maxm);
        ring = ['-sin(psi).*(cos(psi*3)<0.9).*(rad>sin(' num2str(theta_min) '))'];
        [fxc2, fxs2, fyc2, fys2, fzc2, fzs2] = GaussExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, atf, resolution, ring, maxm);
        [fxc2, fxs2, fyc2, fys2, fzc2, fzs2] = RotateEMField(fxc2, fxs2, fyc2, fys2, fzc2, fzs2, -pi/2);
    else
        ring = ['-sin(psi).*(rad>sin(' num2str(theta_min) '))']; maxm = 3;
        [fxc, fxs, fyc, fys, fzc, fzs, rho, z] = GaussExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, atf, resolution, ring, maxm);
        [fxc2, fxs2, fyc2, fys2, fzc2, fzs2] = RotateEMField(fxc, fxs, fyc, fys, fzc, fzs, -pi/2);
    end
end
    
fxc = fxc + fxc2;
fxs = fxs + fxs2;
fyc = fyc + fyc2;
fys = fys + fys2;
fzc = fzc + fzc2;
fzs = fzs + fzs2;

if 0
    subplot(2,1,1)
    FocusImage2D(rho,height-z,cat(3,fxc,fxs))
    colorbar
    subplot(2,1,2)
    FocusImage2D(rho,height-z,cat(3,fzc,fzs))
    colorbar
end

phi = 0:pi/360:2*pi;

if 0
    if eflag % Ez
        subplot(221)
        pcolor(rho(:,1)*cos(phi),rho(:,1)*sin(phi),...
            real((squeeze(fzc(:,1,:))*cos((0:maxm)'*phi) + squeeze(fzs(:,1,:))*sin((1:maxm)'*phi))));
        shading flat; axis image; colorbar('v')
        subplot(222)
        pcolor(rho(:,1)*cos(phi),rho(:,1)*sin(phi),...
            imag((squeeze(fzc(:,1,:))*cos((0:maxm)'*phi) + squeeze(fzs(:,1,:))*sin((1:maxm)'*phi))));
        shading flat; axis image; colorbar('v')
    else % Hphi
        subplot(223)
        pcolor(rho(:,1)*cos(phi),rho(:,1)*sin(phi),...
            real(((squeeze(fxc(:,1,:))*cos((0:maxm)'*phi) + squeeze(fxs(:,1,:))*sin((1:maxm)'*phi)).*(-ones(size(rho,1),1)*sin(phi)) + ...
            (squeeze(fyc(:,1,:))*cos((0:maxm)'*phi) + squeeze(fys(:,1,:))*sin((1:maxm)'*phi)).*(ones(size(rho,1),1)*cos(phi)))));
        shading flat; axis image; colorbar('v')
        subplot(224)
        pcolor(rho(:,1)*cos(phi),rho(:,1)*sin(phi),...
            imag(((squeeze(fxc(:,1,:))*cos((0:maxm)'*phi) + squeeze(fxs(:,1,:))*sin((1:maxm)'*phi)).*(-ones(size(rho,1),1)*sin(phi)) + ...
            (squeeze(fyc(:,1,:))*cos((0:maxm)'*phi) + squeeze(fys(:,1,:))*sin((1:maxm)'*phi)).*(ones(size(rho,1),1)*cos(phi)))));
        shading flat; axis image; colorbar('v')
    end
end

if 0
    subplot(1,3,1)
    pcolor(rho(:,1)*cos(phi),rho(:,1)*sin(phi),real(squeeze(fxc(:,1,:))*cos((0:maxm)'*phi) + squeeze(fxs(:,1,:))*sin((1:maxm)'*phi))); 
    shading flat; axis image
    subplot(1,3,2)
    pcolor(rho(:,1)*cos(phi),rho(:,1)*sin(phi),real(squeeze(fyc(:,1,:))*cos((0:maxm)'*phi) + squeeze(fys(:,1,:))*sin((1:maxm)'*phi))); 
    shading flat; axis image
    subplot(1,3,3)
    pcolor(rho(:,1)*cos(phi),rho(:,1)*sin(phi),real(squeeze(fzc(:,1,:))*cos((0:maxm)'*phi) + squeeze(fzs(:,1,:))*sin((1:maxm)'*phi))); 
    shading flat; axis image
end

if 0
    subplot(1,3,1)
    pcolor(rho(:,1)*cos(phi),rho(:,1)*sin(phi),abs(squeeze(fxc(:,1,:))*cos((0:maxm)'*phi) + squeeze(fxs(:,1,:))*sin((1:maxm)'*phi)).^2); 
    shading flat; axis image
    subplot(1,3,2)
    pcolor(rho(:,1)*cos(phi),rho(:,1)*sin(phi),abs(squeeze(fyc(:,1,:))*cos((0:maxm)'*phi) + squeeze(fys(:,1,:))*sin((1:maxm)'*phi)).^2); 
    shading flat; axis image
    subplot(1,3,3)
    pcolor(rho(:,1)*cos(phi),rho(:,1)*sin(phi),abs(squeeze(fzc(:,1,:))*cos((0:maxm)'*phi) + squeeze(fzs(:,1,:))*sin((1:maxm)'*phi)).^2); 
    shading flat; axis image
end

if 0
    subplot(221)
    pcolor(rho(:,1)*cos(phi),rho(:,1)*sin(phi),...
        real((squeeze(fxc(:,1,:))*cos((0:maxm)'*phi) + squeeze(fxs(:,1,:))*sin((1:maxm)'*phi)).*(ones(size(rho,1),1)*cos(phi)) + ...
        (squeeze(fyc(:,1,:))*cos((0:maxm)'*phi) + squeeze(fys(:,1,:))*sin((1:maxm)'*phi)).*(ones(size(rho,1),1)*sin(phi))));
    shading flat; axis image
    subplot(222)
    pcolor(rho(:,1)*cos(phi),rho(:,1)*sin(phi),...
        imag((squeeze(fxc(:,1,:))*cos((0:maxm)'*phi) + squeeze(fxs(:,1,:))*sin((1:maxm)'*phi)).*(ones(size(rho,1),1)*cos(phi)) + ...
        (squeeze(fyc(:,1,:))*cos((0:maxm)'*phi) + squeeze(fys(:,1,:))*sin((1:maxm)'*phi)).*(ones(size(rho,1),1)*sin(phi))));
    shading flat; axis image
    subplot(223)
    pcolor(rho(:,1)*cos(phi),rho(:,1)*sin(phi),...
        real((squeeze(fxc(:,1,:))*cos((0:maxm)'*phi) + squeeze(fxs(:,1,:))*sin((1:maxm)'*phi)).*(-ones(size(rho,1),1)*sin(phi)) + ...
        (squeeze(fyc(:,1,:))*cos((0:maxm)'*phi) + squeeze(fys(:,1,:))*sin((1:maxm)'*phi)).*(ones(size(rho,1),1)*cos(phi))));
    shading flat; axis image
    subplot(224)
    pcolor(rho(:,1)*cos(phi),rho(:,1)*sin(phi),...
        imag((squeeze(fxc(:,1,:))*cos((0:maxm)'*phi) + squeeze(fxs(:,1,:))*sin((1:maxm)'*phi)).*(-ones(size(rho,1),1)*sin(phi)) + ...
        (squeeze(fyc(:,1,:))*cos((0:maxm)'*phi) + squeeze(fys(:,1,:))*sin((1:maxm)'*phi)).*(ones(size(rho,1),1)*cos(phi))));
    shading flat; axis image
end

if 1
    FocusImage3D(rho,z,cat(4,cat(3,fxc,fxs),cat(3,fyc,fys),cat(3,fzc,fzs)))
    ax = axis;
    axis(ax)
    [fxc, fxs, fyc, fys, fzc, fzs, rho, z] = GaussExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, atf, resolution);
    print -dpng -r300 RadiallyPolarizedExcitationNoObstruction
    figure
    FocusImage3D(rho,z,cat(4,cat(3,fxc,fxs),cat(3,fyc,fys),cat(3,fzc,fzs)))
    axis(ax)
    print -dpng -r300 LinPolForRadiallyPolarizedExcitationNoObstruction
    [fxc2, fxs2, fyc2, fys2, fzc2, fzs2] = RotateEMField(fxc, fxs, fyc, fys, fzc, fzs, pi/2);
    fxc = fxc + 1i*fxc2;
    fxs = fxs + 1i*fxs2;
    fyc = fyc + 1i*fyc2;
    fys = fys + 1i*fys2;
    fzc = fzc + 1i*fzc2;
    fzs = fzs + 1i*fzs2;
    FocusImage3D(rho,z,cat(4,cat(3,fxc,fxs),cat(3,fyc,fys),cat(3,fzc,fzs)))
    axis(ax)
    print -dpng -r300 CircPolForRadiallyPolarizedExcitationNoObstruction
end

if 0
    feldx = FocusImage3D(rho,z,cat(3,fxc,fxs));
    [feldy, psi] = FocusImage3D(rho,z,cat(3,fyc,fys));    
    FocusImage3D(rho,z,abs(feldx.*cos(psi) + feldy.*sin(psi)).^2,1);
end

