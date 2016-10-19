% PSFModel

radv = [1.5 2 3 5]*1e3;
cov = 0;
[x,y] = meshgrid(-2:0.05:2,-2:0.05:2);
rr = sqrt(x.^2+y.^2);
pp = angle(x+i*y);
rhofield = [0 3];
zfield = [-4 4];
NAv = 1.2;
n0 = 1.333;
n = 1.333;
n1 = 1.333;
d0 = [];
d = 0;
d1 = [];
lamex = 0.64;
focpos = 0;
lamem = 0.67;
mag = 60;
avv = [25 50 100];
zpin = 0e3;
atf = [1.52 5];
kappa = 1;
fd = 3e3;
resolution = [20 lamex/0.2];

if 1
    %clear mx my wx wy amp mxl myl wxl wyl ampl
    for jrad = 1:length(radv)
        for jNA = 1:length(NAv);
            for jav = 1:length(avv)
                [jrad jNA jav]
                NA = NAv(jNA);
                over = radv(jrad);
                av = avv(jav);

                exc = GaussExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, atf, resolution);
                mdf = GaussExc2MDF(exc, NA, n0, n, n1, focpos, lamem, mag, av, zpin, atf, kappa);
                rho = exc.rho;
                z = exc.z;
                maxm = (size(mdf.volx,3)-1)/2;  
                
                for jz=1:size(z,2)
                    f = interp1(rho(:,1),mdf.volx(:,jz,1)+mdf.voly(:,jz,1),rr,'cubic');
                    for j=1:maxm
                        f = f + interp1(rho(:,1),mdf.volx(:,jz,j+1)+mdf.voly(:,jz,j+1),rr,'cubic').*cos(j*pp) + ...
                            interp1(rho(:,1),mdf.volx(:,jz,maxm+1+j)+mdf.voly(:,jz,maxm+1+j),rr,'cubic').*sin(j*pp);
                    end
                    [mx(jz,jrad,jNA,jav) ,my(jz,jrad,jNA,jav) ,wx(jz,jrad,jNA,jav) ,wy(jz,jrad,jNA,jav) ,amp(jz,jrad,jNA,jav)] = Gauss2D(x(1,:),y(:,1),f+0.01*max(f(:))*rand(size(f)));

                    fx = interp1(rho(:,1),exc.fxc(:,jz,1),rr,'cubic');
                    fy = interp1(rho(:,1),exc.fyc(:,jz,1),rr,'cubic');
                    fz = interp1(rho(:,1),exc.fzc(:,jz,1),rr,'cubic');
                    for k=1:maxm
                        fx = fx + interp1(rho(:,1),exc.fxc(:,jz,k+1),rr,'cubic').*cos(k*pp) + interp1(rho(:,1),exc.fxs(:,jz,k),rr,'cubic').*sin(k*pp);
                        fy = fy + interp1(rho(:,1),exc.fyc(:,jz,k+1),rr,'cubic').*cos(k*pp) + interp1(rho(:,1),exc.fys(:,jz,k),rr,'cubic').*sin(k*pp);
                        fz = fz + interp1(rho(:,1),exc.fzc(:,jz,k+1),rr,'cubic').*cos(k*pp) + interp1(rho(:,1),exc.fzs(:,jz,k),rr,'cubic').*sin(k*pp);
                    end
                    ff = abs(fx).^2+abs(fy).^2+abs(fz).^2;
                    [mlx(jz,jrad,jNA,jav) ,mly(jz,jrad,jNA,jav) ,wlx(jz,jrad,jNA,jav) ,wly(jz,jrad,jNA,jav) ,ampl(jz,jrad,jNA,jav)] = Gauss2D(x(1,:),y(:,1),ff+0.01*max(ff(:))*rand(size(ff)));
                end
                vpsf(jrad,jNA,jav) = DetectionVolume(exc.rho,exc.z,mdf.volx,mdf.voly);
                ind=abs(z(1,:))<=4;
                las1(:,jrad,jNA,jav) = simplex('LaserBeamFit',[0.4 0],[0 0],[inf 0],[],[],lamex/n,z(1,ind)',sqrt((wx(ind,jrad,jNA,jav).^2+wy(ind,jrad,jNA,jav).^2)/2));
                las2(:,jrad,jNA,jav) = simplex('LaserBeamFit',[0.4 0],[0 0],[inf 0],[],[],lamex/n,z(1,ind)',sqrt((wx(ind,jrad,jNA,jav).^2+wy(ind,jrad,jNA,jav).^2)/2),1,[],amp(ind,jrad,jNA,jav));
                las(:,jrad,jNA,jav) = simplex('LaserBeamFit',[0.4 0],[0 0],[inf 0],[],[],lamex/n,z(1,ind)',sqrt((wx(ind,jrad,jNA,jav).^2+wy(ind,jrad,jNA,jav).^2)/2),1);
                [err, c(jrad,jNA,jav)] = LaserBeamFit(las1(:,jrad,jNA,jav),lamex/n,z(1,ind)',sqrt((wx(ind,jrad,jNA,jav).^2+wy(ind,jrad,jNA,jav).^2)/2));
                cef(:,jrad,jNA,jav) = simplex('D1CEF',[0.1 0],[0 0],[inf 0],[],[],av/mag,lamem/n,z(1,ind)',amp(ind,jrad,jNA,jav));
                vmod(jrad,jNA,jav) = GaussDetectionVolume(las(1,jrad,jNA,jav), cef(1,jrad,jNA,jav), avv(jav)/mag, [lamex lamem]/n);
                eval(['save PSFModelA' mint2str(jrad,2) '-' mint2str(jNA,2) '-' mint2str(jav,2) ' exc mdf NA av jrad jNA jav'])
            end
        end
    end
	eval(['save PSFModelA wx wy amp wlx wly ampl c radv NAv avv'])
end

