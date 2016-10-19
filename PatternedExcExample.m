%close all

lamex = 0.514; % excitation wavelength (in mum)
pixelsize = 0.5; % DMD pixel size in sample plane (in mum)
hadamard_ind = [10 10];% Hadamard index of pattern

nx = 2^7; ny = 2^7; % number of pixels in x- and y-direction
tmp = zeros(nx,ny); 
tmp(hadamard_ind(1),hadamard_ind(2)) = 1; 
tmp = WAT2D(tmp,'hadamard');
tmp = tmp-min(tmp(:));
subdivision = 2^round(log2(pixelsize/lamex*25)); % refinement to achieve ca. lamex/25 grid spacing

pattern = zeros(subdivision*nx, subdivision*ny);
for j=1:size(pattern,1)
    for k=1:size(pattern,2)
        pattern(j,k) = tmp(ceil(j/subdivision),ceil(k/subdivision));
    end
end
[x,y] = meshgrid(-size(pattern,2)/2+0.5:size(pattern,2)/2, -size(pattern,1)/2+0.5:size(pattern,1)/2);
x = x*pixelsize/subdivision;
y = y*pixelsize/subdivision;

NA = 1.2;
fd = 3e3;
n0 = 1.33;
n = 1.33;
n1 = 1.33;
d0 = [];
d = 0; 
d1 = []; 
incoherent = 1;
circular = 0;

zv = 0:0.1:2;

% stepping through excitation intensity distributions along the z-axis
surf_image = 0; % whether to print a surface 3d-plot or a density plot
for jz=1:length(zv)
    zfield = zv(jz);
    [fx, fy, fz] = PatternedExc(x, y, pattern, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, incoherent, circular);
    if surf_image
        if incoherent==0
            surf(x(1:2^9,1:2^9)-min(x(:)),y(1:2^9,1:2^9)-min(y(:)),abs(fx(1:2^9,1:2^9)).^2+abs(fy(1:2^9,1:2^9)).^2);
        else
            surf(x(1:2^9,1:2^9)-min(x(:)),y(1:2^9,1:2^9)-min(y(:)),fx(1:2^9,1:2^9)+fy(1:2^9,1:2^9));
        end
        shading interp; axis tight; colormap jet; view([-20 80])
        title(['\itz\rm = ' int2str(1e3*zfield) ' nm'])
        if jz==1
            ax = axis;
        else
            axis(ax);
        end
    else
        if incoherent==0
            mim(abs(fx(:,:,1)).^2+abs(fy(:,:,1)).^2);
        else
            mim(fx(:,:,1)+fy(:,:,1));
        end
        if jz==1
            ax = caxis;
        else
            caxis(ax);
        end
        title(['\itz\rm = ' int2str(1e3*zfield) ' nm'])
    end
    eval(['print -dpng -r300 Hadamard' mint2str(hadamard_ind(1),3) '_' mint2str(hadamard_ind(2),3) '_z' mint2str(1e3*zfield,4) 'nm']);
end



% everything below this line is not executed and serves only as repository
% for useful code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return 

if incoherent==0
    mim(cat(3,pattern,abs(fx(:,:,1)).^2+abs(fy(:,:,1)).^2));
    caxis([0 max(abs(fx(:)).^2+abs(fy(:).^2))])
else
    mim(cat(3,pattern,fx(:,:,1)+fy(:,:,1)));
    caxis([0 max(abs(fx(:)+fy(:)))])
end

return

z = exc.z(1,:);
subplot(121)
pcolor(x(1,:)'*ones(1,2*length(z)-1),ones(size(x,2),1)*[-fliplr(z(2:end)) z],[fliplr(squeeze(abs(tx(:,round(end/2),2:end).^2))) squeeze(abs(tx(:,round(end/2),:).^2))]/max(abs(tx(:).^2)))
axis image
shading interp
colorbar('v')
xlabel('\itx\rm (\mum)')
ylabel('\itz\rm (\mum)')
subplot(122)
if incoherent==0
    pcolor(x(1,:)'*ones(1,2*length(z)-1),ones(size(x,2),1)*[-fliplr(z(2:end)) z],[fliplr(squeeze(abs(fx(:,round(end/2),2:end).^2))) squeeze(abs(fx(:,round(end/2),:).^2))]/max(abs(fx(:).^2)))
else
    pcolor(x(1,:)'*ones(1,2*length(z)-1),ones(size(x,2),1)*[-fliplr(z(2:end)) z],[fliplr(squeeze(fx(:,round(end/2),2:end))) squeeze(fx(:,round(end/2),:))]/max(fx(:)))
end
axis image
shading interp
colorbar('v')
xlabel('\itx\rm (\mum)')
ylabel('\itz\rm (\mum)')


return

z = exc.z(1,:);
surf(x(1,:)'*ones(1,2*length(z)-1),ones(size(x,2),1)*[-fliplr(z(2:end)) z],[fliplr(squeeze(abs(tx(:,round(end/2),2:end).^2))) squeeze(abs(tx(:,round(end/2),:).^2))]/max(abs(tx(:).^2)))
axis image
shading interp
lighting gouraud

return

z = exc.z(1,:);
pcolor(x(1,:)'*ones(1,2*length(z)-1),ones(size(x,2),1)*[-fliplr(z(2:end)) z],[fliplr(squeeze(abs(fx(:,round(end/2),2:end).^2))) squeeze(abs(fx(:,round(end/2),:).^2))]/max(abs(fx(:).^2)))
axis image
shading interp
colorbar('v')
xlabel('\itx\rm (\mum)')
ylabel('\itz\rm (\mum)')

return

z = exc.z(1,:);
subplot(121)
pcolor(x(1,:)'*ones(1,2*length(z)-1),ones(size(x,2),1)*[-fliplr(z(2:end)) z],[fliplr(squeeze(abs(tx(:,round(end/2),2:end).^2))) squeeze(abs(tx(:,round(end/2),:).^2))]/max(abs(tx(:).^2)))
axis image
shading interp
colorbar('v')
xlabel('\itx\rm (\mum)')
ylabel('\itz\rm (\mum)')
subplot(122)
pcolor(x(1,:)'*ones(1,2*length(z)-1),ones(size(x,2),1)*[-fliplr(z(2:end)) z],[fliplr(squeeze(abs(fx(:,round(end/2),2:end).^2))) squeeze(abs(fx(:,round(end/2),:).^2))]/max(abs(fx(:).^2)))
axis image
shading interp
colorbar('v')
xlabel('\itx\rm (\mum)')
ylabel('\itz\rm (\mum)')

return 

tsh = 1./exp(1:3);
tsh = sort(tsh);
if length(tsh)>1
    al = fliplr(1./(1:length(tsh)));
    co = flipud([ones(length(tsh),1) (0:length(tsh)-1)'/length(tsh) zeros(length(tsh),1)]);
else
    al = 0.5;
    co = [1 0.8 0.5];
end

z = exc.z(1,:);
tmp = abs(fx.^2)+abs(permute(fx,[2 1 3]).^2);
close all
for j=2:length(z)
    for k=1:length(tsh)
        patch(isosurface(repmat(x,[1 1 2]),repmat(y,[1 1 2]),permute(repmat(z(j-1:j)',[1, size(x,1) size(x,2)]),[2 3 1]),tmp(:,:,j-1:j),tsh(k)*max(tmp(:))),'FaceColor',co(k,:),'EdgeColor','none','FaceAlpha',al(k));
    end
    if j==2
        axis image
        view(3);
        material metal
        light('position',[0 1 0])
        light('position',[1 -1 0])
        light('position',[-1 -1 0])
        lighting gouraud
        axis vis3d
        box on
        axis([min(x(:)) max(x(:)) min(y(:)) max(y(:)) min(z) max(z)])
        cameratoolbar
    end
    %drawnow
    eval(['print -dpng -r300 Checkerboard' mint2str(j-1,2)])
end
