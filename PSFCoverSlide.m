% PSFModel

close all

if 1
    radv = (1:0.05:1.8)*1e3;
    covv = -5:5;
    abev = 0:0.01:0.05;
    [x,y] = meshgrid(-2:0.05:2,-2:0.05:2);
    rr = sqrt(x.^2+y.^2);
    pp = angle(x+i*y);
    rhofield = [0 3];
    zfield = [-4 4];
    NA = 1.14;
    n0 = 1.333;
    n = 1;
    n1 = 1;
    d0 = [];
    d = 0;
    d1 = [];
    lamex = 0.64;

    clear mx1 my1 wx1 wy1 amp1 mx2 my2 wx2 wy2 amp2
    % load PSFModel_NA114
    for jrad = 1:length(radv)
        for jcov = 1:length(covv)
            for jabe = 1:length(abev)
                [jrad jcov jabe]
                abe = abev(jabe);
                over = [0 radv(jrad) 1];
                dfoc = 0.205;
                focpos = [[0 dfoc 0 0 0];[0 -dfoc 0 0 0]];
                pow = 1;
                lamem = 0.67;
                mag = 60;
                av = 100;
                zpin = 0e3;
                atf = [1.52 covv(jcov)];
                kappa = 1;,
                lt = [];
                sat = 0;
                fd = 3e3;
                resolution = [20 lamex/0.2];
                ring = [];
                maxm = 10;

                [volx1, voly1, volx2, voly2, rho, z, fxc1, fxs1, fyc1, fys1, fzc1, fzs1, fxc2, fxs2, fyc2, fys2, fzc2, fzs2] = DICPSF(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, pow, lamem, mag, av, zpin, atf, kappa, lt, sat, resolution, ring, abe, maxm);

                eval(['save DICPSF' mint2str(jrad,2) mint2str(jcov,2) mint2str(jabe,2) ' volx1 voly1 volx2 voly2 rho z fxc1 fxs1 fyc1 fys1 fzc1 fzs1 fxc2 fxs2 fyc2 fys2 fzc2 fzs2'])

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

                    [mx1(jz,jrad,jcov,jabe),my1(jz,jrad,jcov,jabe),wx1(jz,jrad,jcov,jabe),wy1(jz,jrad,jcov,jabe),amp1(jz,jrad,jcov,jabe)] = Gauss2D(x(1,:),y(:,1),f1+0.01*max(f1(:))*rand(size(f1)));
                    [mx2(jz,jrad,jcov,jabe),my2(jz,jrad,jcov,jabe),wx2(jz,jrad,jcov,jabe),wy2(jz,jrad,jcov,jabe),amp2(jz,jrad,jcov,jabe)] = Gauss2D(x(1,:),y(:,1),f2+0.01*max(f2(:))*rand(size(f1)));
                end
            end
        end
    end
    eval(['save PSFModel_NA114'])

    return

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

        [mx1a(jz,kk),my1a(jz,kk),wx1a(jz,kk),wy1a(jz,kk),amp1a(jz,kk)] = Gauss2D(x(1,:),y(:,1),f1+0.1*max(f1(:)));
        [mx2a(jz,kk),my2a(jz,kk),wx2a(jz,kk),wy2a(jz,kk),amp2a(jz,kk)] = Gauss2D(x(1,:),y(:,1),f2+0.1*max(f2(:)));
    end


    return
    
    for jrad = 1:length(radv);
        for jast = 1:length(astv);
            wx = wx1(maxl(-wx1(:,jrad)),jrad);
            wy = wy1(maxl(-wy1(:,jrad)),jrad);            
            if length(wx)==1 & length(wy)==1
                w0x1(jrad) = wx;
                w0y1(jrad) = wy;                
                z0x1(jrad) = z(1,wx1(:,jrad)==wx);
                z0y1(jrad) = z(1,wy1(:,jrad)==wy);
            elseif length(wx)==2
                w0x1(jrad) = wx(1);
                w0y1(jrad) = wx(2);                
                z0x1(jrad) = z(1,wx1(:,jrad)==wx(1));
                z0y1(jrad) = z(1,wx1(:,jrad)==wx(2));
            else
                w0x1(jrad) = wy(1);
                w0y1(jrad) = wy(2);                
                z0x1(jrad) = z(1,wy1(:,jrad)==wy(1));
                z0y1(jrad) = z(1,wy1(:,jrad)==wy(2));
            end
            wx = wx2(maxl(-wx2(:,jrad)),jrad);
            wy = wy2(maxl(-wy2(:,jrad)),jrad);            
            if length(wx)==1 & length(wy)==1
                w0x2(jrad) = wx;
                w0y2(jrad) = wy;                
                z0x2(jrad) = z(1,wx2(:,jrad)==wx);
                z0y2(jrad) = z(1,wy2(:,jrad)==wy);
            elseif length(wx)==2
                w0x2(jrad) = wx(1);
                w0y2(jrad) = wx(2);                
                z0x2(jrad) = z(1,wx2(:,jrad)==wx(1));
                z0y2(jrad) = z(1,wx2(:,jrad)==wx(2));
            else
                w0x2(jrad) = wy(1);
                w0y2(jrad) = wy(2);                
                z0x2(jrad) = z(1,wy2(:,jrad)==wy(1));
                z0y2(jrad) = z(1,wy2(:,jrad)==wy(2));
            end
        end
    end
    
    return
    
    load D:\Joerg\Doc\Fcs\2Focus\PSF\PSF_2005-12-06\PSFFitResults beamwaist beampos
    
    close all
    for jast=1:length(astv)
        figure
        plot(radv/1e3,1e3*w0x1(:),radv/1e3,1e3*w0y1(:)); ax = axis;
        plot(radv/1e3,1e3*w0x1(:),radv/1e3,1e3*w0y1(:),ax(1:2),beamwaist(1,1)*[1 1],':',ax(1:2),beamwaist(2,1)*[1 1],':');
        title(['Foc 1 Astig = ' num2str(astv(jast))]);
        figure
        plot(radv/1e3,1e3*w0x2(:),radv/1e3,1e3*w0y2(:)); ax = axis;
        plot(radv/1e3,1e3*w0x2(:),radv/1e3,1e3*w0y2(:),ax(1:2),beamwaist(3,1)*[1 1],':',ax(1:2),beamwaist(4,1)*[1 1],':');
        title(['Foc 2 Astig = ' num2str(astv(jast))]);
    end
    
    return

    load D:\Joerg\Doc\Fcs\2Focus\PSF\PSF_2005-12-06\PSFResults.mat
    ex.wx1=wx1; ex.wx2=wx2; ex.wy1=wy1; ex.wy2=wy2; ex.mx1=mx1;ex.mx2=mx2;ex.my1=my1;ex.my2=my2;ex.d=d;ex.z=zd;
    ex.om1 = om1; ex.om2 = om2; ex.amp1=amp1; ex.amp2=amp2;

    tst = simplex('LaserBeamFit',[500 100],[],[],[],[],640/1.33,zd,wy1.*(zd(1,:)>0)+wx1.*(zd(1,:)<=0));
    [err,c] = LaserBeamFit(tst,640/1.33,zd,wy1.*(zd(1,:)>0)+wx1.*(zd(1,:)<=0));
    bw = tst(1)*c;
    tst = simplex('LaserBeamFit',[500 100],[],[],[],[],640/1.33,zd,wx1.*(zd(1,:)>0)+wy1.*(zd(1,:)<=0));
    [err,c] = LaserBeamFit(tst,640/1.33,zd,wx1.*(zd(1,:)>0)+wy1.*(zd(1,:)<=0));
    bw = [bw tst(1)*c];
    tst = simplex('LaserBeamFit',[500 100],[],[],[],[],640/1.33,zd,wy2.*(zd(1,:)>=0)+wx2.*(zd(1,:)<0));
    [err,c] = LaserBeamFit(tst,640/1.33,zd,wy2.*(zd(1,:)>=0)+wx2.*(zd(1,:)<0));
    bw = [bw tst(1)*c];
    tst = simplex('LaserBeamFit',[500 100],[],[],[],[],640/1.33,zd,wx2.*(zd(1,:)>=0)+wy2.*(zd(1,:)<0));
    [err,c] = LaserBeamFit(tst,640/1.33,zd,wx2.*(zd(1,:)>=0)+wy2.*(zd(1,:)<0));
    bw = [bw tst(1)*c];

    close all
    for jast=1:length(astv)
        figure
        plot(radv/1e3,1e3*w0x1(:),radv/1e3,1e3*w0y1(:)); ax = axis;
        plot(radv/1e3,1e3*w0x1(:),radv/1e3,1e3*w0y1(:),ax(1:2),bw(1)*[1 1],':',ax(1:2),bw(2)*[1 1],':');
        title(['Foc 1 Astig = ' num2str(astv(jast))]);
        figure
        plot(radv/1e3,1e3*w0x2(:),radv/1e3,1e3*w0y2(:)); ax = axis;
        plot(radv/1e3,1e3*w0x2(:),radv/1e3,1e3*w0y2(:),ax(1:2),bw(3)*[1 1],':',ax(1:2),bw(4)*[1 1],':');
        title(['Foc 2 Astig = ' num2str(astv(jast))]);
    end
    
    return
    
    load D:\Joerg\Doc\Fcs\2Focus\PSF\PSF_2006-02-20\PSFexp.mat
    
    
    plot(ex.z-250,ex.wx2,'o',z(1,:)*1e3,wx1(:,6,2)*1e3,ex.z-250,ex.wy1,'o',z(1,:)*1e3,wy2(:,6,2)*1e3)
    plot(ex.z-200,ex.wx2,'o',-z(1,:)*1e3,wx2(:,8,3)*1e3,ex.z-200,ex.wy2,'o',-z(1,:)*1e3,wy1(:,8,3)*1e3)
    
    x1 = interp1(-z(1,:)*1e3,1e3*mx1(:,8,3),ex.z-150,'cubic');
    y1 = interp1(-z(1,:)*1e3,1e3*my1(:,8,3),ex.z-150,'cubic');
    x2 = interp1(-z(1,:)*1e3,1e3*mx2(:,8,3),ex.z-200,'cubic');
    y2 = interp1(-z(1,:)*1e3,1e3*my2(:,8,3),ex.z-200,'cubic');
    d = sqrt((x2-x1).^2+(y2-y1).^2);
    plot(ex.z,ex.d,ex.z,d)
    
end


if 0
    radv = (1:0.1:2)*1e3;
    astv = 0:0.1:0.5;
    cov = 0;
    [x,y] = meshgrid(-2:0.05:2,-2:0.05:2);
    rr = sqrt(x.^2+y.^2);
    pp = angle(x+i*y);
    rhofield = [0 3];
    zfield = [-4 4];
    NA = 1.14;
    n0 = 1.333;
    n = n0;
    n1 = n0;
    d0 = [];
    d = 0;
    d1 = [];
    lamex = 0.64;

    clear mx my wx wy om amp
    for jrad = 1:length(radv);
        for jast = 1:length(astv);
            over = [astv(jast) radv(jrad) 1];
            dfoc = 0.205;
            focpos = 0;
            lamem = 0.67;
            mag = 60;
            av = 100;
            zpin = 0e3;
            atf = [];
            kappa = 0;,
            lt = [];
            sat = 0;
            fd = 3e3;
            resolution = [20 lamex/0.2];
            ring = [];
            maxm = 10;

            [volx, voly, rho, z] = OneFocusMDF(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, lamem, mag, av, zpin, atf, kappa, lt, sat, resolution, ring, maxm);

            for jz=1:size(z,2)
                f1 = interp1(rho(:,1),volx(:,jz,1)+voly(:,jz,1),rr,'cubic');
                for j=1:maxm
                    f1 = f1 + interp1(rho(:,1),volx(:,jz,j+1)+voly(:,jz,j+1),rr,'cubic').*cos(j*pp) + ...
                        interp1(rho(:,1),volx(:,jz,maxm+1+j)+voly(:,jz,maxm+1+j),rr,'cubic').*sin(j*pp);
                end

                [mx(jz,jrad),my(jz,jrad),wx(jz,jrad),wy(jz,jrad),om(jz,jrad),amp(jz,jrad)] = GaussEllipse(x(1,:),y(:,1),f1+0.1*max(f1(:)));
            end
        end
    end
    eval(['save OneFocusPSF_NA114'])
end


if 0
    rad = (1:0.1:2)*1e3;
    ast = -0.3:0.01:0.3;
    cov = -5:5;
    [x,y] = meshgrid(-2:0.05:2,-2:0.05:2);
    rr = sqrt(x.^2+y.^2);
    pp = angle(x+i*y);
    rhofield = [0 2];
    zfield = [-4 4];
    NA = 1.14;
    n0 = 1.333;
    n = 1.333;
    n1 = 1.333;
    d0 = [];
    d = 0;
    d1 = [];
    lamex = 0.635;

    mx1 = zeros(size(z,2),length(ast),length(cov)); my1 = mx1; wx1 = mx1; wy1 = mx1; om1 = mx1; amp1 = mx1;
    mx2 = mx1; my2 = mx1; wx2 = mx1; wy2 = mx1; om2 = mx1; amp2 = mx1;
    for jrad = 1:length(rad);
        for jast = 1:length(ast)
            for jcov = 1:length(cov)
                over = [ast(jast) rad(jrad) 1];
                dfoc = 0.2/sqrt(2);
                focpos = [0 dfoc dfoc 0 0];
                pow = 1;
                lamem = 0.67;
                mag = 60;
                av = 75;
                zpin = 0e3;
                atf = [1.52 cov(jcov)];
                kappa = 1;,
                lt = [];
                sat = 0;
                fd = 3e3;
                resolution = [20 lamex/0.2];
                ring = [];
                maxm = 10;

                [volx1, voly1, volx2, voly2, rho, z] = TwoFocusMDF(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, pow, lamem, mag, av, zpin, atf, kappa, lt, sat, resolution, ring, maxm);

                for jz=1:size(z,2)
                    f1 = interp1(rho(:,1),volx1(:,jz,1)+voly1(:,jz,1),rr,'cubic');
                    for j=1:maxm
                        f1 = f1 + interp1(rho(:,1),volx1(:,jz,j+1)+voly1(:,jz,j+1),rr,'cubic').*cos(j*pp) + ...
                            interp1(rho(:,1),volx1(:,jz,maxm+1+j)+voly1(:,jz,maxm+1+j),rr,'cubic').*sin(j*pp);
                    end
                    f2 = interp1(rho(:,1),volx1(:,jz,1)+voly1(:,jz,1),rr,'cubic');
                    for j=1:maxm
                        f2 = f2 + interp1(rho(:,1),volx2(:,jz,j+1)+voly2(:,jz,j+1),rr,'cubic').*cos(j*pp) + ...
                            interp1(rho(:,1),volx2(:,jz,maxm+1+j)+voly2(:,jz,maxm+1+j),rr,'cubic').*sin(j*pp);
                    end

                    [mx1(jz,jcov),my1(jz,jcov),wx1(jz,jcov),wy1(jz,jcov),om1(jz,jcov),amp1(jz,jcov)] = GaussEllipse(x(1,:),y(:,1),f1+0.1*max(f1(:)));
                    [mx2(jz,jcov),my2(jz,jcov),wx2(jz,jcov),wy2(jz,jcov),om2(jz,jcov),amp2(jz,jcov)] = GaussEllipse(x(1,:),y(:,1),f2+0.1*max(f2(:)));
                end
            end
        end
        eval(['save PSFModel_rad' mint2str(over(2),3)])
    end
end

