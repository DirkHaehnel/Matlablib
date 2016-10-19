function [orient, imm, err] = RefineOrientation(im,thetav,phiv,nn,pixel,rho,fxx0,fxx2,fxz,byx0,byx2,byz);

bck = Disk(nn);
bck = bck/sqrt(sum(sum(bck.^2)));

for j=1:length(thetav)
    al = thetav(j);
    for k=1:length(phiv)
        be = -phiv(k);
        mask = SEPImage(al,be,nn,pixel,rho,fxx0,fxx2,fxz,byx0,byx2,byz).*(bck>0);
        mask = mask./sqrt(sum(sum(mask.^2)));
        crs = sum(sum(bck.*mask));
        crs = inv([1 crs; crs 1]);
        for kx=1:(size(im,2)-2*nn)
            for ky=1:(size(im,1)-2*nn)
                tmp = im(ky:ky+2*nn,kx:kx+2*nn);
                im0 = sum(sum(tmp.*bck));
                im02 = sum(sum(tmp.^2.*bck/max(bck(:))));
                im1 = sum(sum(tmp.*mask));
                tmperr(j,k,kx,ky) = im02 - crs(1,1)*im0*im0 - 2*crs(1,2)*im0*im1 - crs(2,2)*im1*im1;
            end
        end
    end
end
phii = phiv(1):diff(phiv([1 end]))/40:phiv(end);
thetai = thetav(1):diff(thetav([1 end]))/40:thetav(end);
[x,y] = meshgrid(1:length(phii),1:length(thetai));
for kx=1:(size(im,2)-2*nn)
    for ky=1:(size(im,1)-2*nn)
        erri = interp2(phiv,thetav,tmperr(:,:,kx,ky),phii,thetai','cubic');
        ind1 = x(erri==min(erri(:)));
        ind2 = y(erri==min(erri(:)));
        orient(2,kx,ky) = phii(ind1(1));
        orient(1,kx,ky) = thetai(ind2(1));
        imm(:,:,kx,ky) = bck.*SEPImage(orient(1,kx,ky),-orient(2,kx,ky),nn,pixel,rho,fxx0,fxx2,fxz,byx0,byx2,byz);
        err(kx,ky) = erri(ind2(1),ind1(1));
    end
end