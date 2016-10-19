% program for modeling the imaging of inclined tubules in TIRFM under
% non-physical conditions ignoring all interfaceeffects in imaging

clear all
close all

lamex = 0.491;
lamem = 0.525;
n0 = 1.52;
n = n0;
n1 = n0;
d0 = [];
d = 10;
d1 = [];
NA = 1.46;
mag = 100;
fd = 164.5e3/mag;
pixel = 0.1;
av = sqrt(pixel^2/pi)*mag;
resolution = [lamex/0.01 lamex/0.02]; % spatial grid resolution for MDF calculation 
focpos = 0; % position of focus in mum
rhofield =  [0 3*d];
zfield = focpos + [-lamex/resolution(2)/2 d];
over = [0 0];
atf = [];
zpin = 0;
kappa = 0;
lt = 1;

% calculating the detection efficiency distribution
exc.NA = NA;
exc.fd = fd;
exc.n0 = n0;
exc.n = n;
exc.n1 = n1;
exc.d0 = d0/2/pi*lamex;
exc.d = d/2/pi*lamex;
exc.d1 = d1/2/pi*lamex;
exc.focpos = focpos;
exc.atf = atf;
exc.ring = [];
exc.maxm = 0;

drho = lamex/resolution(1);
dz = lamex/resolution(2);
rhov = rhofield(1) + (0.5:(rhofield(2)-rhofield(1))/drho)*drho;
zv = zfield(1) + (0.5:(zfield(2)-zfield(1))/dz)*dz;
[exc.rho, exc.z] = ndgrid(rhov, zv);
exc.fxc = ones(size(exc.rho));
exc.fyc = 0*exc.fxc;
exc.fzc = 0*exc.fxc;
exc.fxs = [];
exc.fys = 0*exc.fxs;
exc.fzs = 0*exc.fxs;

wv = 0.05:0.05:2;
incv = (3:12)/180*pi;

mdf = GaussExc2MDF(exc, NA, n0, n, n1, focpos, lamem, mag, av, zpin, atf, kappa, 1);
% FocusImage3D(exc.rho,exc.z,mdf.volx+mdf.voly);


rho = exc.rho(:,1)';
psf = zeros(size(rho,1),size(exc.z,2));
for j=1:length(exc.z(1,:))
    for k=1:length(rho)
        psf(k,j) = rho(k:end)./sqrt(rho(k:end).^2-(rho(k)-rho(1))^2)*(mdf.volx(k:end,j,1) + mdf.voly(k:end,j,1));
    end
    disp(j)
end

for jw = 1:length(wv) 
    psf1 = psf.*exp(-exc.z(1:length(rho),:)/wv(jw));
    for jinc = 1:length(incv)
        tic
        % tubulus image
        inclinationangle = incv(jinc);
        len = 3*wv(jw)/sin(inclinationangle);
        xt = drho/2:drho:len*cos(inclinationangle); 
        zt = xt*tan(inclinationangle);
        xi = 0:pixel:len;
        im = 0*xi; im1 = im;
        for j=1:length(xt)
            zz = zt(j)/dz;
            ind = ceil(zz);
            tmp = psf1(:,ind)*(ind-zz) + (zz-ind+1)*psf1(:,ind+1);
            im1 = im1 + interp1(rho,tmp,abs(xi-xt(j)),'cubic',0);
            tmp = psf(:,ind)*(ind-zz) + (zz-ind+1)*psf(:,ind+1);
            im = im + interp1(rho,tmp,abs(xi-xt(j)),'cubic',0);
        end

        ind = 5:length(im)-5;
        c = polyfit(xi(1,ind)*tan(inclinationangle),log(im1(ind)./im(ind)),1);
        p1(jinc,jw) = -1/c(1)/wv(jw);

        disp([jw jinc toc])
    end
end

save InclinedTubulusInTIRFIdeal

wi = wv(1)+ (0:200)/200*diff(wv([1 end]));
inci = pi/180*(3:0.05:12);
surf(wi, inci/pi*180, interp2(ones(length(incv),1)*wv,incv'*ones(1,length(wv)),p1,ones(length(inci),1)*wi,inci'*ones(1,length(wi)),'cubic'))
shading interp
cameratoolbar
xlabel('penetration depth \itd\rm (\mum)'); ylabel('inclination angle (°)'); zlabel('\itd_{fit}\rm/\itd\rm')


