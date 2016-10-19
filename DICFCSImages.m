radv = [1.25 1.5 2 3 4]*1e3;
covv = 0:10;
astv = 0:0.1:0.3;

dist = 0.40;
rhofield = [0 1.5];
zfield = [9 23];
NA = 1.14;
fd = 3e3;
n0 = 1.333;
n = n0;
n1 = n0;
d0 = [];
d = 0;
d1 = [];
lamex = 0.64;
focposexc = [15 dist/2 0 0 0];
pow = 1;
lamem = 0.67;
mag = 60;
av = 100;
focposdet = 15;
zpin = 0e3;
kappa = 0;
lt = [];
pulse = [0.05 25]/2; % laser characteristics in units of lifetime
triplet = 0;
resolution = [30 10];
ring = [];
maxm = 10;

satv = [0 0.05*exp(log(1/0.05)*(0:9)/9)];


for jrad = 3:3
    for jcov = 1:1
        for jsat = 11:11

            over = radv(jrad);
            atf = [1.52 covv(jcov)];
            tic
            exc = DICExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focposexc, pow, atf, resolution, maxm);
            mdf = DICExc2MDF(exc, NA, n0, n, n1, focposdet, lamem, mag, av, zpin, atf, kappa, lt, pulse, satv(jsat), triplet);

            phi=0:pi/100:2*pi;

            mm = max(squeeze(mdf.volx1(1,:,1)+mdf.voly1(1,:,1)));
            tst = squeeze(mdf.volx1(1,:,1)+mdf.voly1(1,:,1))/mm;
            z = exc.z(1,:)-15;
            rho = exc.rho(:,1);
            ind = 1:length(z);
            % ind = [min(ind(tst>0.1)) min(ind(tst>0.3)) min(ind(tst>0.8)) max(ind(tst==1)) max(ind(tst>0.8)) max(ind(tst>0.3)) max(ind(tst>0.1))];
            % ind = [max(ind(tst==1)) max(ind(tst>0.8)) max(ind(tst>0.3)) max(ind(tst>0.1))];
            % if length(ind)>4 ind(1)=[]; end
            ind = [max(ind(z<0))+1 max(ind(tst>0.1))];
            ind = round(ind(1)+diff(ind)*[0 1/3 2/3 1]);

            FocusImage3D(exc.rho,exc.z-15,mdf.volx1+mdf.voly1);
            hold on
            for j=1:length(ind)
                surf(exc.rho(:,1)*cos(phi),exc.rho(:,1)*sin(phi),z(ind(j))*ones(size(exc.rho,1),length(phi)),...
                    permute(repmat([0 0 1]',[1 size(exc.rho,1) length(phi)]),[2 3 1]),'facealpha',0.2,'edgecolor','none')
            end
            axis image
            hold off

            eval(['print -dpng -r300 DICFCSMDF' mint2str(jrad,2) 'rad_' mint2str(jcov,2) 'cov_' mint2str(jsat,2) 'sat'])

            clear ex fx p y
            for j=1:length(ind)
                fx(:,j) = FocusImage2D(rho,z,mdf.volx1(:,ind(j),:)+mdf.volx1(:,ind(j),:),0,'calc');
                ex(:,j) = FocusImage2D(rho,z,mdf.volx1(:,ind(j),:)+mdf.volx1(:,ind(j),:),pi,'calc');
                p(1,j) = [-flipud(rho);rho]'*[flipud(ex(:,j));fx(:,j)]./sum([flipud(ex(:,j));fx(:,j)]);
                p(2,j) = sqrt([-flipud(rho);rho].^2'*[flipud(ex(:,j));fx(:,j)]./sum([flipud(ex(:,j));fx(:,j)])-p(1,j).^2)/2;
                for k=1:2
                    p(:,j) = Simplex('Gauss',p(:,j),[],[],[],[],[-rho(end:-1:1);rho],[ex(end:-1:1,j); fx(:,j)],0,[],0,[ex(end:-1:1,j); fx(:,j)]);
                end
                [err, c, y(:,j)] = Gauss(p(:,j),[-rho(end:-1:1);rho],[ex(end:-1:1,j); fx(:,j)],0);
            end
            plot([-flipud(rho);rho],[flipud(ex);fx]/mm,'o',[-flipud(rho);rho],y/mm)
            colorize
            xlabel('\itx\rm [\mum]'); ylabel('norm. intensity')
            for j=1:length(ind)
                s{j} = [mnum2str(z(ind(j)),1,1) ' \mum'];
            end
            legend(s,2)
            eval(['print -dpng -r300 DICFCSGauss' mint2str(jrad,2) 'rad_' mint2str(jcov,2) 'cov_' mint2str(jsat,2) 'sat'])

        end
    end
end
