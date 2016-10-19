function [feld, phi, rr, zz] = FocusImage3D(rho,z,vol,flag,tsh,maxangle,clfflag,varargin)

if nargin<5 || isempty(tsh)
    tsh = 1./exp(1:3);
end
if nargin<6 || isempty(maxangle)
    maxangle = 2*pi;
end

if ndims(rho)>2
    rr = rho;
    zz = z;
    feld = vol;
    phi = flag;
else
    maxphi = 2^7;
    [zz,rr,phi] = meshgrid(z(1,:),rho(:,1),(0:maxphi)/maxphi*maxangle);

    if nargin<4 || isempty(flag)
        if nargout>0
            feld = zeros(size(rr,1),size(rr,2),size(rr,3),size(vol,4));
        else
            feld = 0*rr;
        end
        for k=1:size(rr,3)
            psi = squeeze(phi(1,1,k));
            tmp = vol(:,:,1,:);
            for j=1:(size(vol,3)-1)/2
                tmp = tmp + vol(:,:,j+1,:)*cos(j*psi) + vol(:,:,(end+1)/2+j,:)*sin(j*psi);
            end
            if nargout==0
                if isreal(vol)
                    feld(:,:,k) = sum(tmp,4);
                else
                    feld(:,:,k) = sum(abs(tmp).^2,4);
                end
            else
                feld(:,:,k,:) = tmp;
            end
        end
    else
        feld = sum(vol,4);
    end
end

if nargout==0
    tsh = sort(tsh);
    if length(tsh)>1
        al = fliplr(1./(1:length(tsh)));
        co = flipud([ones(length(tsh),1) (0:length(tsh)-1)'/length(tsh) zeros(length(tsh),1)]);
    else
        al = 0.5;
        co = [1 0.8 0.5];
    end

    if nargin<7 || isempty(clfflag)
        clf;
    end
    for j=1:length(tsh)
        patch(isosurface(rr.*cos(phi),rr.*sin(phi),zz,feld,tsh(j)*max(feld(:))),'FaceColor',co(j,:),'EdgeColor','none','FaceAlpha',al(j));
    end

    axis image
    view(3);
    material metal
    if nargin<7 || isempty(clfflag)
        light('position',[0 1 0])
        light('position',[1 -1 0])
        light('position',[-1 -1 0])
    end
    lighting gouraud
    axis vis3d
    box on
    if nargin>7
        set(gca,varargin{:});
    end
    cameratoolbar
    
    xlabel('\itx\rm (\mum)'); ylabel('\ity\rm (\mum)'); zlabel('\itz\rm (\mum)')

    clear feld phi
end

% mpcolor(squeeze(rr(:,50,:).*cos(phi(:,50,:))),squeeze(rr(:,50,:).*sin(phi(:,50,:))),squeeze(feld(:,50,:)))