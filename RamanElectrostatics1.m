% program Raman electrostatics

close all

probe_height = 50;
probe_protrusion_distance = 5;
probe_radius = 5;
protrusion_radius = 5;
width = 50;
lateral_shift = 5;
height = protrusion_radius + probe_protrusion_distance + probe_height;
delta = 1;

x = (-width + min(0, lateral_shift)):delta:(width + max(0, lateral_shift));
y = -width:delta:width;
z = 0:delta:height;

[x, y, z] = meshgrid(x, y, z);

a = sqrt(probe_radius*probe_height);

probe = sqrt(((x-lateral_shift).^2 + y.^2)/a^2 + (z - height).^2/probe_height^2) <= 1;
protrusion = sqrt(x.^2+y.^2+z.^2) <= protrusion_radius;

psi = z/height;
psi(protrusion) = 0;
psi(probe) = 1;

ly = size(psi,1);
lx = size(psi,2);
lz = size(psi,3);
len = numel(psi);
e = ones(len,1);
M = spdiags([e/6 e/6 e/6 -e e/6 e/6 e/6],[-ly*lx -ly -1 0 1 ly ly*lx],len,len);
boundary = false(numel(psi),1);
boundary(x(:)==min(x(:))) = true;
boundary(x(:)==max(x(:))) = true;
boundary(y(:)==min(y(:))) = true;
boundary(y(:)==max(y(:))) = true;
boundary(z(:)==min(z(:))) = true;
boundary(z(:)==max(z(:))) = true;
boundary(probe(:)) = true;
boundary(protrusion(:)) = true;

psi(~boundary) = M(~boundary,~boundary)\(M(~boundary,boundary)*psi(boundary));

contourf(squeeze(x((end+1)/2,:,:)),squeeze(z((end+1)/2,:,:)),squeeze(psi((end+1)/2,:,:)),0:0.025:1);
axis image
drawnow

[ex, ey, ez] = gradient(psi);
figure
pcolor(squeeze(x((end+1)/2,:,:)),...
    squeeze(z((end+1)/2,:,:)),...
    log10(squeeze(ex((end+1)/2,:,:).^2 + ey((end+1)/2,:,:).^2 + ez((end+1)/2,:,:).^2)/delta^2*height^2));
axis image
colorbar
drawnow

return


