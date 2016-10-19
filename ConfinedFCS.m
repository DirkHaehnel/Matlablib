% program for studying FCS in a spherical confinement

radiusv = [0.5 1 2 3];
rv = (0:60)/60;
pv = 0:pi/36:pi/2;

resolution = 60;
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

modres = FCS(exc.rho,exc.z,mdf.volx,mdf.voly); modres=[1; modres];

%clear res
load D:\Joerg\MATLAB\ConfinedFCSRes.mat
for jr=4:length(radiusv)
    radius = radiusv(jr); % cell radius in mum

    Nsub = 3;
    lam = 2^(1/Nsub);
    timewin = 1e7;
    Ncasc = ceil(log2(timewin));
    autotime = [0;lam.^(0:Ncasc*Nsub-1)'];
    dif = 5e-5; % diffusion constant in mum^2/time unit

    maxzeros = 50;
    clear kalpha
    for j=1:maxm+1
        kalpha(j,:) = mBesselZero(j-1,1e-20,maxzeros,3)/radius;
        if j>1 && kalpha(j,1)<1e-5
            tmp = mBesselZero(j-1,1e-20,maxzeros+1,3)/radius;
            kalpha(j,:) = tmp(2:end);
        end
    end

    for kr=1:length(rv)
        for jp=1:length(pv)
            if ~(kr==1 && jp>1)
                x0 = (radius+2)*rv(kr)*sin(pv(jp));
                z0 = (radius+2)*rv(kr)*cos(pv(jp));
                maxm = max([10 ceil(sqrt(x0.^2+z0.^2)*25)]);

                maxtheta = min([pi max([max(max(angle(z0-exc.z+i*(x0+exc.rho)))) max(max(angle(z0-exc.z+i*(x0-exc.rho))))])]);
                if maxtheta<=0
                    maxtheta = pi;
                end
                mintheta = max([0 min([min(min(angle(z0+exc.z+i*(x0-exc.rho)))) min(min(angle(z0+exc.z+i*(x0+exc.rho))))])]);
                if x0<=max(exc.rho(:,1))
                    maxphi = pi;
                else
                    maxphi = asin(max(exc.rho(:,1))/x0);
                end
                minr = min([0 min(min(sqrt((x0-exc.rho).^2+(z0-exc.z).^2)))]);
                maxr = min([radius max(max(sqrt((x0+exc.rho).^2+(z0+exc.z).^2)))]);

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
                    clear nrm
                    for j=1:maxm+1
                        for k=1:maxzeros
                            if j==1 && k==1
                                nrm(j,k) = 2/3*radius^2;
                            else
                                nrm(j,k) = pi*besselj(j-0.5,kalpha(j,k)*radius)^2/kalpha(j,k)/radius/2*(radius^2-j*(j-1)/kalpha(j,k)^2);
                            end
                        end
                    end
                    clear Jv
                    for j=1:maxm+1
                        for k=1:maxzeros
                            if j==1 && k==1
                                Jv(:,k,j) = rad.^2.*sqrt(2/pi);
                            else
                                Jv(:,k,j) = sqrt(1/kalpha(j,k).*rad.^3).*besselj(j-0.5,kalpha(j,k)*rad);
                            end
                        end
                    end

                    auto = 0*autotime;
                    coef = zeros(maxm+1,maxzeros,maxm+1);
                    intensity(kr,jp,jr) = 0;
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
                                coef(j,:,kk+1) = coef(j,:,kk+1) + (2-(jj==0))/sqrt(1+(kk==0))*((Jv(:,:,j)'*tmp*M{j}(kk+1,:)')*cos(kk*psi))';
                            end
                        end
                        jj
                        toc
                    end
                    
                    intensity(kr,jp,jr) = coef(1,1,1)/(2*maxm+1);
                    coef = sum(coef.^2,3)./nrm;
                    for j=1:maxm+1
                        for k=1:maxzeros
                            auto = auto + coef(j,k)*exp(-dif*kalpha(j,k).^2*autotime);
                        end
                    end
                    auto = auto/auto(1);
                    res(:,kr,jp,jr) = auto;
                else
                    res(:,kr,jp,jr) = 0;
                end
            end
        end
    end
    for jp=2:length(pv)
        res(:,1,jp,jr) = res(:,1,1,jr);
    end
    save ConfinedFCSRes res exc mdf modres autotime radiusv pv rv intensity
    
    for kr=1:size(res,2)
        for jp=1:size(res,3)
            for jr=1:size(res,4)
                if sum(res(:,kr,jp,jr))>0 && all(isfinite(res(:,kr,jp,jr)))
                    tau(kr,jp,jr) = interp1(res(diff(res(:,kr,jp,jr))<0,kr,jp,jr),autotime(diff(res(:,kr,jp,jr))<0),0.5*(1+res(end,kr,jp,jr)),'cubic');
                else
                    tau(kr,jp,jr) = 0;
                end
            end
        end
    end
    tau0 = interp1(modres(diff(modres)<0),autotime(diff(modres)<0),0.5*(1+modres(end)),'cubic');
    save ConfinedFCSRes tau tau0 -append

end

return

j=3; 
tau(1,2:end,j)=tau(1,1,j);
pcolor((2+radiusv(j))*rv'*sin(pv),(2+radiusv(j))*rv'*cos(pv),tau(:,:,j)/tau0); hold on; pcolor(-(2+radiusv(j))*rv'*sin(pv),(2+radiusv(j))*rv'*cos(pv),tau(:,:,j)/tau0); pcolor((2+radiusv(j))*rv'*sin(pv),-(2+radiusv(j))*rv'*cos(pv),tau(:,:,j)/tau0); pcolor(-(2+radiusv(j))*rv'*sin(pv),-(2+radiusv(j))*rv'*cos(pv),tau(:,:,j)/tau0); axis image;
shading interp
plot(radiusv(j)*sin(0:pi/100:2*pi),radiusv(j)*cos(0:pi/100:2*pi),'y'); hold off
colorbar
patch([-0.5 0.5 0.5 -0.5],[-radiusv(j)-1.55 -radiusv(j)-1.55 -radiusv(j)-1.45 -radiusv(j)-1.45],'w')
axis off
%set(gca,'xticklabel',abs(get(gca,'xtick')))

