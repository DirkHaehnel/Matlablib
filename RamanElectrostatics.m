% program Raman electrostatics

close all

probe_height = 50;
probe_protrusion_distance = 5;
probe_radius = 5;
protrusion_radius = 5;
width = 50;
lateral_shift = 5;
height = protrusion_radius + probe_protrusion_distance + probe_height;
delta = 0.25;

x = (-width + min(0, lateral_shift)):delta:(width + max(0, lateral_shift));
y = -width:delta:width;
z = 0:delta:height;

[x, y, z] = meshgrid(x, y, z);

a = sqrt(probe_radius*probe_height);

probe = sqrt(((x-lateral_shift).^2 + y.^2)/a^2 + (z - height).^2/probe_height^2) <= 1;
protrusion = sqrt(x.^2+y.^2+z.^2) <= protrusion_radius;

load metals
eps_au = 1.5^2; %gold(wavelength==514).^2;

psi = z/height;
psi(protrusion) = 0;
psi(probe) = 1;
% eps = ones(size(psi));
% eps(protrusion) = eps_au;
% eps(probe) = eps_au;
% laplace_eps = (eps([1 1:end-1],:,:) + eps([2:end end],:,:) + ...
%         eps(:,[1 1:end-1],:) + eps(:,[2:end end],:) + ...
%         eps(:,:,[1 1:end-1]) + eps(:,:,[2:end end]))/6 - eps;

gauss = 0.8;
tst = (psi([1 1:end-1],:,:) + psi([2:end end],:,:) + ...
        psi(:,[1 1:end-1],:) + psi(:,[2:end end],:) + ...
        psi(:,:,[1 1:end-1]) + psi(:,:,[2:end end]))/6 - psi;
while sum(abs(tst(:)))>1e-10
    %epspsi = eps.*psi;
    psi = psi + gauss*tst;
    %     psi = psi + gauss*(eps.*((psi([1 1:end-1],:,:) + psi([2:end end],:,:) + ...
    %         psi(:,[1 1:end-1],:) + psi(:,[2:end end],:) + ...
    %         psi(:,:,[1 1:end-1]) + psi(:,:,[2:end end]))/6 - psi) + ...
    %         (epspsi([1 1:end-1],:,:) + epspsi([2:end end],:,:) + ...
    %         epspsi(:,[1 1:end-1],:) + epspsi(:,[2:end end],:) + ...
    %         epspsi(:,:,[1 1:end-1]) + epspsi(:,:,[2:end end]))/6 - epspsi - ...
    %         psi.*laplace_eps);
    psi(protrusion) = 0;
    psi(probe) = 1;
    psi(:,:,end) = 1;
    psi(:,:,1) = 0;    
    psi(1,:,:) = z(1,:,:)/height;
    psi(end,:,:) = z(end,:,:)/height;
    psi(:,1,:) = z(:,1,:)/height;
    psi(:,end,:) = z(:,end,:)/height;
    tst = (psi([1 1:end-1],:,:) + psi([2:end end],:,:) + ...
        psi(:,[1 1:end-1],:) + psi(:,[2:end end],:) + ...
        psi(:,:,[1 1:end-1]) + psi(:,:,[2:end end]))/6 - psi;

    %pcolor(squeeze(x((end+1)/2,:,:)),squeeze(z((end+1)/2,:,:)),squeeze(psi((end+1)/2,:,:)))
    %shading interp
    contourf(squeeze(x((end+1)/2,:,:)),squeeze(z((end+1)/2,:,:)),squeeze(psi((end+1)/2,:,:)),0:0.025:1); 
    axis image
    drawnow
end

[ex, ey, ez] = gradient(psi);
pcolor(squeeze(x((end+1)/2,:,:)),...
    squeeze(z((end+1)/2,:,:)),...
    log10(squeeze(ex((end+1)/2,:,:).^2 + ey((end+1)/2,:,:).^2 + ez((end+1)/2,:,:).^2)/delta^2*height^2));
axis image
colorbar
drawnow

return


