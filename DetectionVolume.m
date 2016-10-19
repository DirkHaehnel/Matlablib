function [veff, intens] = DetectionVolume(rho,z,volx,voly)

if nargin<4 || isempty(voly)
   voly = volx; 
end
nn = size(volx,3);
rhom = repmat(rho,[1 1 nn]);
dz = z(1,2)-z(1,1); drho = rho(2,1)-rho(1,1); 

volx = volx/max(volx(:));
voly = voly/max(voly(:));

intens = zeros(1,1+size(voly,4));
intens(1) = squeeze(sum(sum(rhom(:,:,1).*volx(:,:,1))))*drho*dz*pi*2;
for j=1:size(voly,4)
    intens(j+1) = squeeze(sum(sum(rhom(:,:,1).*voly(:,:,1,j))))*drho*dz*pi*2;
end

tmp = voly(:,:,:,1);
veff = zeros(size(voly,4),1);
for j=1:size(voly,4)
    veff(j) = real([2 ones(1,nn-1)]*squeeze(sum(sum(rhom.*volx.*tmp))))*(drho*dz*pi);
    if j<size(voly,4)
        tmp = FieldConv(tmp,voly(:,:,:,j+1));
    end
end

