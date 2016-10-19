% PSFModel

radv = (1:0.1:4)*1e3;
cov = 0;
abev = [];
[x,y] = meshgrid(-2:0.05:2,-2:0.05:2);
rr = sqrt(x.^2+y.^2);
pp = angle(x+i*y);
rhofield = [0 3];
zfield = [-4 4];
NAv = 1:0.1:1.3;
n0 = 1.333;
n = 1.333;
n1 = 1.333;
d0 = [];
d = 0;
d1 = [];
lamex = 0.64;
dfoc = 0.205;
focpos = [[0 dfoc 0 0 0];[0 -dfoc 0 0 0]];
pow = 1;
lamem = 0.67;
mag = 60;
av = 100;
zpin = 0e3;
atf = [];
kappa = 1;,
lt = [];
sat = 0;
fd = 3e3;
resolution = [20 lamex/0.2];
ring = [];
maxm = 10;

if 0
    clear mx1 my1 wx1 wy1 amp1 mx2 my2 wx2 wy2 amp2
    for jrad = 1:length(radv)
        for jNA = 1:length(NAv);
            [jrad jNA]
            NA = NAv(jNA);
            over = [0 radv(jrad) 1];

            [volx1, voly1, volx2, voly2, rho, z, fxc1, fxs1, fyc1, fys1, fzc1, fzs1, fxc2, fxs2, fyc2, fys2, fzc2, fzs2] = DICMDF(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, pow, lamem, mag, av, zpin, atf, kappa, lt, sat, resolution, ring, maxm);

            % eval(['save DICPSF' mint2str(jrad,2) mint2str(jcov,2) mint2str(jabe,2) ' volx1 voly1 volx2 voly2 rho z fxc1 fxs1 fyc1 fys1 fzc1 fzs1 fxc2 fxs2 fyc2 fys2 fzc2 fzs2'])

            for jz=1:size(z,2)
                f1 = interp1(rho(:,1),volx1(:,jz,1)+voly1(:,jz,1),rr,'cubic');
                for j=1:maxm
                    f1 = f1 + interp1(rho(:,1),volx1(:,jz,j+1)+voly1(:,jz,j+1),rr,'cubic').*cos(j*pp) + ...
                        interp1(rho(:,1),volx1(:,jz,maxm+1+j)+voly1(:,jz,maxm+1+j),rr,'cubic').*sin(j*pp);
                end
                f2 = interp1(rho(:,1),volx2(:,jz,1)+voly2(:,jz,1),rr,'cubic');
                for j=1:maxm
                    f2 = f2 + interp1(rho(:,1),volx2(:,jz,j+1)+voly2(:,jz,j+1),rr,'cubic').*cos(j*pp) + ...
                        interp1(rho(:,1),volx2(:,jz,maxm+1+j)+voly2(:,jz,maxm+1+j),rr,'cubic').*sin(j*pp);
                end

                [mx1(jz,jrad,jNA),my1(jz,jrad,jNA),wx1(jz,jrad,jNA),wy1(jz,jrad,jNA),amp1(jz,jrad,jNA)] = Gauss2D(x(1,:),y(:,1),f1+0.01*max(f1(:))*rand(size(f1)));
                [mx2(jz,jrad,jNA),my2(jz,jrad,jNA),wx2(jz,jrad,jNA),wy2(jz,jrad,jNA),amp2(jz,jrad,jNA)] = Gauss2D(x(1,:),y(:,1),f2+0.01*max(f2(:))*rand(size(f1)));


                fx = interp1(rho(:,1),fxc1(:,jz,1),rr,'cubic');
                fy = interp1(rho(:,1),fyc1(:,jz,1),rr,'cubic');
                fz = interp1(rho(:,1),fzc1(:,jz,1),rr,'cubic');
                for k=1:maxm
                    fx = fx + interp1(rho(:,1),fxc1(:,jz,k+1),rr,'cubic').*cos(k*pp) + interp1(rho(:,1),fxs1(:,jz,k),rr,'cubic').*sin(k*pp);
                    fy = fy + interp1(rho(:,1),fyc1(:,jz,k+1),rr,'cubic').*cos(k*pp) + interp1(rho(:,1),fys1(:,jz,k),rr,'cubic').*sin(k*pp);
                    fz = fz + interp1(rho(:,1),fzc1(:,jz,k+1),rr,'cubic').*cos(k*pp) + interp1(rho(:,1),fzs1(:,jz,k),rr,'cubic').*sin(k*pp);
                end
                f1 = abs(fx).^2+abs(fy).^2+abs(fz).^2;

                fx = interp1(rho(:,1),fxc2(:,jz,1),rr,'cubic');
                fy = interp1(rho(:,1),fyc2(:,jz,1),rr,'cubic');
                fz = interp1(rho(:,1),fzc2(:,jz,1),rr,'cubic');
                for k=1:maxm
                    fx = fx + interp1(rho(:,1),fxc2(:,jz,k+1),rr,'cubic').*cos(k*pp) + interp1(rho(:,1),fxs2(:,jz,k),rr,'cubic').*sin(k*pp);
                    fy = fy + interp1(rho(:,1),fyc2(:,jz,k+1),rr,'cubic').*cos(k*pp) + interp1(rho(:,1),fys2(:,jz,k),rr,'cubic').*sin(k*pp);
                    fz = fz + interp1(rho(:,1),fzc2(:,jz,k+1),rr,'cubic').*cos(k*pp) + interp1(rho(:,1),fzs2(:,jz,k),rr,'cubic').*sin(k*pp);
                end
                f2 = abs(fx).^2+abs(fy).^2+abs(fz).^2;

                [mlx1(jz,jrad,jNA),mly1(jz,jrad,jNA),wlx1(jz,jrad,jNA),wly1(jz,jrad,jNA),ampl1(jz,jrad,jNA)] = Gauss2D(x(1,:),y(:,1),f1+0.01*max(f1(:))*rand(size(f1)));
                [mlx2(jz,jrad,jNA),mly2(jz,jrad,jNA),wlx2(jz,jrad,jNA),wly2(jz,jrad,jNA),ampl2(jz,jrad,jNA)] = Gauss2D(x(1,:),y(:,1),f2+0.01*max(f2(:))*rand(size(f1)));
            end
        end
    end
    eval(['save PSFModel'])
end

if 1
    %load PSFModel
    ind=abs(z(1,:))<2;
    for j=1:size(wx1,2)
        for k=1:size(wx1,3)
            las1(:,j,k) = simplex('LaserBeamFit',[0.4 0],[],[],[],[],lamex/n,z(1,ind)',sqrt((wx1(ind,j,k).^2+wy1(ind,j,k).^2)/2));
            [err, c1(j,k)] = LaserBeamFit(las1(:,j,k),lamex/n,z(1,ind)',sqrt((wx1(ind,j,k).^2+wy1(ind,j,k).^2)/2));
            las2(:,j,k) = simplex('LaserBeamFit',[0.4 0],[],[],[],[],lamex/n,z(1,ind)',sqrt((wx2(ind,j,k).^2+wy2(ind,j,k).^2)/2));
            [err, c2(j,k)] = LaserBeamFit(las2(:,j,k),lamex/n,z(1,ind)',sqrt((wx2(ind,j,k).^2+wy2(ind,j,k).^2)/2));
            cef1(:,j,k) = simplex('D1CEF',[1 0 0.2],[],[],[],[],av/mag,lamem/n,z(1,ind)',amp1(ind,j,k)./wx1(ind,j,k)./wy1(ind,j,k));
            cef2(:,j,k) = simplex('D1CEF',[1 0 0.2],[],[],[],[],av/mag,lamem/n,z(1,ind)',amp2(ind,j,k)./wx2(ind,j,k)./wy2(ind,j,k));
        end
    end
end


