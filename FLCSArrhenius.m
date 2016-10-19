function [err, c, z] = FLCSArrehnius(p, t, y, Temp, bld);

% Function [err, c, z] = FLCSArrehnius(p, t, y) calculates the error err, amplitudes c, and fit curve z for the 
% FCS data y with the model offset + 1/(p(1)+t)/sqrt(p(2)+t)*(1 + exp(-t/p(3)) + exp(-t/p(4)) + exp(-t/p(5)))
% input y has to be a n-by-4 matrix

if length(t)<size(y,1)
    c = t;
    t = y;
    p = p(:)';
    a = p(1);
    b = p(2);
    om = reshape(p(3:2+4*(size(y,4)-1)),4,size(y,4)-1); om = [ones(4,1) om];
    p = p(3+4*(size(y,4)-1):end);
    fac = p(1:3)';
    arh = p(4:6)';
    fac = (fac*ones(1,size(y,4))).*exp(-arh*(1./Temp(:)'));
    col = ones(size(t));
    
    for k=1:size(y,4)
        for j=1:4
            zz = [ones(length(t),1) ((a*sqrt(b)*om(j,k)./(a*om(j,k)+t)./sqrt(b*om(j,k)+t))*ones(1,3)).*[col exp(-t*fac(1:2,k)')]];
            zm = (a*sqrt(b)*om(j,k)./(a*om(j,k)+t)./sqrt(b*om(j,k)+t)).*exp(-t*fac(3,k));
            switch j
                case 1
                    c(:,1,1,k) = lsqnonneg([zz zm],y(:,1,1,k));
                    z(:,1,1,k) = [zz zm]*c(:,1,1,k);
                case 2
                    c(:,1,2,k) = lsqnonneg([zz -zm],y(:,1,2,k)); c(end,1,2,k) = -c(end,1,2,k);
                    z(:,1,2,k) = [zz zm]*c(:,1,2,k);
                case 3
                    c(:,2,1,k) = lsqnonneg([zz -zm],y(:,2,1,k)); c(end,2,1,k) = -c(end,2,1,k);
                    z(:,2,1,k) = [zz zm]*c(:,2,1,k);
                case 4
                    c(:,2,2,k) = lsqnonneg([zz zm],y(:,2,2,k));
                    z(:,2,2,k) = [zz zm]*c(:,2,2,k);
            end
        end
    end
    
    err = z;
else
    t = t(:);
    p = p(:)';
    a = p(1);
    b = p(2);
    om = reshape(p(3:2+4*(size(y,4)-1)),4,size(y,4)-1); om = [ones(4,1) om];
    p = p(3+4*(size(y,4)-1):end);
    fac = p(1:3)';
    arh = p(4:6)';
    fac = (fac*ones(1,size(y,4))).*exp(-arh*(1./Temp(:)'));
    col = ones(size(t));
    
    for k=1:size(y,4)
        for j=1:4
            zz = [ones(length(t),1) ((a*sqrt(b)*om(j,k)./(a*om(j,k)+t)./sqrt(b*om(j,k)+t))*ones(1,3)).*[col exp(-t*fac(1:2,k)')]];
            zm = (a*sqrt(b)*om(j,k)./(a*om(j,k)+t)./sqrt(b*om(j,k)+t)).*exp(-t*fac(3,k));
            switch j
                case 1
                    c(:,1,1,k) = lsqnonneg([zz zm],y(:,1,1,k));
                    z(:,1,1,k) = [zz zm]*c(:,1,1,k);
                case 2
                    c(:,1,2,k) = lsqnonneg([zz -zm],y(:,1,2,k)); c(end,1,2,k) = -c(end,1,2,k);
                    z(:,1,2,k) = [zz zm]*c(:,1,2,k);
                case 3
                    c(:,2,1,k) = lsqnonneg([zz -zm],y(:,2,1,k)); c(end,2,1,k) = -c(end,2,1,k);
                    z(:,2,1,k) = [zz zm]*c(:,2,1,k);
                case 4
                    c(:,2,2,k) = lsqnonneg([zz zm],y(:,2,2,k));
                    z(:,2,2,k) = [zz zm]*c(:,2,2,k);
            end
        end
    end
    clear zz
    for k=1:size(y,4)
        for j1=1:2
            for j2=1:2
                yy(:,4*(k-1)+(j2-1)*2+j1) = y(:,j1,j2,k);
                zz(:,4*(k-1)+(j2-1)*2+j1) = z(:,j1,j2,k);
            end
        end
    end
    if nargin>4 & ~isempty(bld)
        semilogx(t,yy./(col*mean(yy)),'--',t,zz./(col*mean(zz))); drawnow;
    end
    err = sum(sum(yy-zz,2).^2./abs(sum(zz,2)));
end

