% program for studying FCS in a spherical confinement

radiusv = [0.1 0.2 0.5 1];
rmax = 10;
rv = (0:60)/60;
pv = 0:pi/36:pi/2;

resolution = 100;
lamex = 0.635;

rhofield = [-lamex/resolution/2 0.5+lamex/resolution];
zfield = [-lamex/resolution/2 1+lamex/resolution];
NA = 1.2;
fd = 5e3;
n0 = 1.333;
n = 1.333;
n1 = 1.333;
d0 = [];
d = [];
d1 = [];
over = 5e3;
focpos = 0;
atf = [];
ring = '';
lamem = 0.67;
mag = 60;
av = 10; %pinhole radius!!!
zpin = 0;

exc = GaussExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, atf, resolution);
mdf = GaussExc2MDF(exc, NA, n0, n, n1, focpos, lamem, mag, av, zpin, atf);

clear res
for jr=1:1 %1:length(radiusv)
    radius = radiusv(jr); % cell radius in mum

    Nsub = 3;
    lam = 2^(1/Nsub);
    timewin = 1e7;
    Ncasc = ceil(log2(timewin));
    autotime = [0;lam.^(0:Ncasc*Nsub-1)'];
    dif = 5e-5; % diffusion constant in mum^2/time unit

    kalpha = 0.5:0.5:30;%10.^(-1:0.1:log10(20));
    maxzeros = length(kalpha);

    for kr=30:30 %length(rv)
        for jp=3:3%19:19 %length(pv)
            if 1 %~(kr==1 && jp>1)
                x0 = (radius+2)*rv(kr)*sin(pv(jp));
                z0 = (radius+2)*rv(kr)*cos(pv(jp));
                maxm = max([10 ceil(sqrt(x0.^2+z0.^2)*35)]);

                tmp = [z0-exc.z+i*(x0+exc.rho) z0-exc.z+i*(x0-exc.rho)];
                maxtheta = min([pi max(max((angle(tmp).*(abs(tmp)>=radius).*(imag(tmp)>=0))))]);
                %maxtheta = min([pi max([max(max(angle(z0-exc.z+i*(x0+exc.rho)))) max(max(angle(z0-exc.z+i*(x0-exc.rho))))])]);
                if maxtheta<=0
                    maxtheta = pi;
                end
                tmp = [z0+exc.z+i*(x0+exc.rho) z0+exc.z+i*(x0-exc.rho)];
                mintheta = max([0 min(min((angle(tmp).*(abs(tmp)>=radius).*(imag(tmp)>=0))))]);
                %mintheta = max([0 min([min(min(angle(z0+exc.z+i*(x0-exc.rho)))) min(min(angle(z0+exc.z+i*(x0+exc.rho))))])]);
                if x0<=max(exc.rho(:,1))
                    maxphi = pi;
                else
                    maxphi = asin(max(exc.rho(:,1))/x0);
                end
                minr = max([radius min(min(sqrt((x0-exc.rho).^2+(z0-exc.z).^2)))]);
                maxr = max(max(sqrt((x0+exc.rho).^2+(z0+exc.z).^2)));

                rad = (minr:lamex/resolution:maxr)';
                if ~isempty(rad)
                    da = pi/360;
                    theta = round(mintheta/da)*da:da:round(maxtheta/da)*da;
                    [tt,rr] = meshgrid(theta,rad);
                    cc = cos(tt); ss = sin(tt);
                    clear M
                    for j=1:maxm+1
                        M{j} = legendre(j-1,cc(1,:),'norm').*(ones(j,1)*ss(1,:));
                    end
                    clear Jv
                    for j=1:maxm+1
                        for k=1:maxzeros
                            Jv(:,k,j) = sqrt(1/kalpha(k).*rad.^3).*besselj(j-0.5,kalpha(k)*rad);
                            tmp = sqrt(1/kalpha(k).*rad.^3).*bessely(j-0.5,kalpha(k)*rad);
                            tstj = sqrt(1/kalpha(k)/radius).*(-besselj(j-0.5,kalpha(k)*radius)./radius + ...
                                kalpha(k)*(besselj(j-1.5,kalpha(k)*radius)-besselj(j+0.5,kalpha(k)*radius)));
                            tsty = sqrt(1/kalpha(k)/radius).*(-bessely(j-0.5,kalpha(k)*radius)./radius + ...
                                kalpha(k)*(bessely(j-1.5,kalpha(k)*radius)-bessely(j+0.5,kalpha(k)*radius)));
                            %Jv(:,k,j) = tsty*Jv(:,k,j) - tstj*tmp;
                            Jv(:,k,j) = (tsty*Jv(:,k,j) - tstj*tmp)/sqrt(tsty^2+tstj^2)*kalpha(k);
                        end
%                         for k=1:maxzeros
%                             for kk=1:k-1
%                                 Jv(:,k,j) = Jv(:,k,j) - Jv(:,kk,j)*sum(Jv(:,k,j).*Jv(:,kk,j)./rad.^2);
%                             end
%                             Jv(:,k,j) = Jv(:,k,j)/sqrt(sum(Jv(:,k,j).^2./rad.^2));
%                         end
                    end

                    coef = zeros(maxm+1,maxzeros,maxm+1);
                    for jj=0:maxm % mirror symmetry around xz-plane
                        tic
                        psi = jj/(2*maxm+1)*2*pi;
                        rho = sqrt((rr.*ss*cos(psi)-x0).^2 + (rr.*ss*sin(psi)).^2);
                        z = rr.*cc-z0;
                        ind = rho<=max(exc.rho(:,1)) & abs(z)<max(exc.z(1,:));
                        tmp = zeros(length(rad),length(theta));
                        if sum(ind(:))>0
                            tmp(ind) = griddata(exc.rho,exc.z,mdf.volx(:,:,1)+mdf.voly(:,:,1),rho(ind),abs(z(ind)));
                        end
                        for j=1:maxm+1
                            for kk=0:j-1
                                coef(j,:,kk+1) = coef(j,:,kk+1) + (2-(jj==0))/(1+(kk==0))*((Jv(:,:,j)'*tmp*M{j}(kk+1,:)')*cos(kk*psi))';
                            end
                        end
                        jj
                        toc
                    end

                    bild = 0*tt;
                    for j=1:maxm+1
                        for k=1:maxzeros
                            for kk=0:j-1
                                bild = bild + coef(j,k,kk+1)*cos(kk*0)*(Jv(:,k,j)./(rad.^2+(rad==0)))*(M{j}(kk+1,:)./(ss(1,:)+(ss(1,:)==0)));
                            end
                        end
                    end

                    pcolor(rr.*sin(tt),rr.*cos(tt),bild); axis image; shading interp
                else
                end
            end
        end
    end
end
