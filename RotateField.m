function [fxcr, fxsr, fycr, fysr, fzcr, fzsr] = RotateField(fxc, fxs, fyc, fys, fzc, fzs, phi)

if nargin==2
    if isfield(fxc,'fxc')
        phi = fxs;
        tmp = fxc;
        tmp.fxc = cos(phi)*fxc.fxc-sin(phi)*fxc.fyc;
        tmp.fyc = sin(phi)*fxc.fxc+cos(phi)*fxc.fyc;
        tmp.fzc = fxc.fzc;
        tmp.fxs = fxc.fys;
        tmp.fys = fxc.fys;
        tmp.fzs = fxc.fzs;
        if ~isempty(fxc.fxs)
            for j=1:size(fxc.fxs,3)
                tmp.fxc(:,:,j+1) = -sin(phi)*(fxc.fyc(:,:,j+1)*cos(j*phi) - fxc.fys(:,:,j)*sin(j*phi)) + cos(phi)*(fxc.fxc(:,:,j+1)*cos(j*phi) - fxc.fxs(:,:,j)*sin(j*phi));
                tmp.fxs(:,:,j) = -sin(phi)*(fxc.fyc(:,:,j+1)*sin(j*phi) + fxc.fys(:,:,j)*cos(j*phi)) + cos(phi)*(fxc.fxc(:,:,j+1)*sin(j*phi) + fxc.fxs(:,:,j)*cos(j*phi));
                tmp.fyc(:,:,j+1) = sin(phi)*(fxc.fxc(:,:,j+1)*cos(j*phi) - fxc.fxs(:,:,j)*sin(j*phi)) + cos(phi)*(fxc.fyc(:,:,j+1)*cos(j*phi) - fxc.fys(:,:,j)*sin(j*phi));
                tmp.fys(:,:,j) = sin(phi)*(fxc.fxc(:,:,j+1)*sin(j*phi) + fxc.fxs(:,:,j)*cos(j*phi)) + cos(phi)*(fxc.fyc(:,:,j+1)*sin(j*phi) + fxc.fys(:,:,j)*cos(j*phi));
                tmp.fzc(:,:,j+1) = (fxc.fzc(:,:,j+1)*cos(j*phi) - fxc.fzs(:,:,j)*sin(j*phi));
                tmp.fzs(:,:,j) = (fxc.fzc(:,:,j+1)*sin(j*phi) + fxc.fzs(:,:,j)*cos(j*phi));
            end
        end
        fxcr = tmp;
    elseif isfield(fxc,'volx') && isfield(fxc,'voly')
        phi = fxs;
        tmp = fxc;
        tmp.volx = fxc.volx;
        tmp.voly = fxc.voly;
        for j=1:(size(fxc.volx,3)-1)/2;
            tmp.volx(:,:,j+1) = fxc.volx(:,:,j+1)*cos(j*phi) - fxc.volx(:,:,(end+1)/2+j)*sin(j*phi);
            tmp.volx(:,:,(end+1)/2+j) = fxc.volx(:,:,j+1)*sin(j*phi) + fxc.volx(:,:,(end+1)/2+j)*cos(j*phi);
            tmp.voly(:,:,j+1) = fxc.voly(:,:,j+1)*cos(j*phi) - fxc.voly(:,:,(end+1)/2+j)*sin(j*phi);
            tmp.voly(:,:,(end+1)/2+j) = fxc.voly(:,:,j+1)*sin(j*phi) + fxc.voly(:,:,(end+1)/2+j)*cos(j*phi);
        end
        fxcr = tmp;
    elseif isfield(fxc,'volx') 
        phi = fxs;
        tmp = fxc;
        tmp.volx = fxc.volx;
        for j=1:(size(fxc.volx,3)-1)/2;
            tmp.volx(:,:,j+1) = fxc.volx(:,:,j+1)*cos(j*phi) - fxc.volx(:,:,(end+1)/2+j)*sin(j*phi);
            tmp.volx(:,:,(end+1)/2+j) = fxc.volx(:,:,j+1)*sin(j*phi) + fxc.volx(:,:,(end+1)/2+j)*cos(j*phi);
        end
        fxcr = tmp;
    end        
else
    fxcr = cos(phi)*fxc-sin(phi)*fyc;
    fycr = sin(phi)*fxc+cos(phi)*fyc;
    fzcr = fzc;
    fxsr = fys;
    fysr = fys;
    fzsr = fzs;
    if ~isempty(fxs)
        for j=1:size(fxs,3)
            fxcr(:,:,j+1) = -sin(phi)*(fyc(:,:,j+1)*cos(j*phi) - fys(:,:,j)*sin(j*phi)) + cos(phi)*(fxc(:,:,j+1)*cos(j*phi) - fxs(:,:,j)*sin(j*phi));
            fxsr(:,:,j) = -sin(phi)*(fyc(:,:,j+1)*sin(j*phi) + fys(:,:,j)*cos(j*phi)) + cos(phi)*(fxc(:,:,j+1)*sin(j*phi) + fxs(:,:,j)*cos(j*phi));
            fycr(:,:,j+1) = sin(phi)*(fxc(:,:,j+1)*cos(j*phi) - fxs(:,:,j)*sin(j*phi)) + cos(phi)*(fyc(:,:,j+1)*cos(j*phi) - fys(:,:,j)*sin(j*phi));
            fysr(:,:,j) = sin(phi)*(fxc(:,:,j+1)*sin(j*phi) + fxs(:,:,j)*cos(j*phi)) + cos(phi)*(fyc(:,:,j+1)*sin(j*phi) + fys(:,:,j)*cos(j*phi));
            fzcr(:,:,j+1) = (fzc(:,:,j+1)*cos(j*phi) - fzs(:,:,j)*sin(j*phi));
            fzsr(:,:,j) = (fzc(:,:,j+1)*sin(j*phi) + fzs(:,:,j)*cos(j*phi));
        end
    end
end