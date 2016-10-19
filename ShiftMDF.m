function mdfs = ShiftMDF(mdf, delta, maxm)

mdfs = mdf;
ind = mdf.rho(:,1) + abs(delta) < max(mdf.rho(:,1)); 
mdfs.rho = mdfs.rho(ind,:);
mdfs.z = mdfs.z(ind,:);
mdfs.volx = zeros(sum(ind),size(mdf.z,2),2*maxm+1);
mdfs.voly = mdfs.volx;
rho = mdf.rho(:,1);
row = ones(1,size(mdf.z,2));
psi = (0:4*maxm)'/(4*maxm+1)*2*pi;
mc = cos(psi*(1:maxm));
ms = sin(psi*(1:maxm));
for jr=1:size(mdfs.rho,1)
    phi = unwrap(angle(-delta + rho(jr)*(cos(psi) + 1i*sin(psi))));
    rad = sqrt((-delta + rho(jr)*cos(psi)).^2 + (rho(jr)*sin(psi)).^2);

    tmp = interp1(rho,mdf.volx(:,:,1),rad,'cubic','extrap');
    for j = 1:(size(mdf.volx,3)-1)/2
        tmp = tmp + interp1(rho,mdf.volx(:,:,j+1),rad,'cubic','extrap').*(cos(j*phi)*row) + interp1(rho,mdf.volx(:,:,(end+1)/2+j),rad,'cubic','extrap').*(sin(j*phi)*row);
    end
    if sum(isnan(tmp))>0
        disp('stop')
    end
    mdfs.volx(jr,:,1) = sum(tmp);
    mdfs.volx(jr,:,2:maxm+1) = (2*mc'*tmp)';
    mdfs.volx(jr,:,maxm+2:2*maxm+1) = (2*ms'*tmp)';
    mdfs.volx(jr,:,:) = mdfs.volx(jr,:,:)/(4*maxm+1);
    
    tmp = interp1(rho,mdf.voly(:,:,1),rad,'cubic','extrap');
    for j = 1:(size(mdf.voly,3)-1)/2
        tmp = tmp + interp1(rho,mdf.voly(:,:,j+1),rad,'cubic','extrap').*(cos(j*phi)*row) + interp1(rho,mdf.voly(:,:,(end+1)/2+j),rad,'cubic','extrap').*(sin(j*phi)*row);
    end
    mdfs.voly(jr,:,1) = sum(tmp);
    mdfs.voly(jr,:,2:maxm+1) = (2*mc'*tmp)';
    mdfs.voly(jr,:,maxm+2:2*maxm+1) = (2*ms'*tmp)';
    mdfs.voly(jr,:,:) = mdfs.voly(jr,:,:)/(4*maxm+1);
end

