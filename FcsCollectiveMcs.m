% Program for modeling FCS on interacting particles 

close all

if ~exist('u','var')
    lamex = 0.543; % excitation wavelength
    resolution = [lamex/0.025 lamex/0.05]; % spatioal grid resolution for MDF calculation 
    rhomax = 3; % radial extent of calculation volume in mum
    zmax = 5; % axial extent of calculation volume in mum
    rhofield =  [0 sqrt(2)*rhomax + (lamex/resolution(1))/2];
    focpos = 45; % position of focus in mum
    zfield = focpos + [0 zmax+(lamex/resolution(2))/2];
    NA = 0.9; % numerical aperture
    n0 = 1.497; % ref. index of solution/immersion medium
    n = n0;
    n1 = n0;
    d0 = [];
    d = 0;
    d1 = [];
    mag = 40; % magnification
    fd = 164.5e3/mag; % focal distance of objective in mum
    over = 5e3; % radius of laser beam in mum
    ring = [];
    av = 70/2; % radius of confocal aperture in mum
    lamem = 0.570; % center emission wavelength
    zpin = 0; % position of confocal aperture

    exc = GaussExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, [], resolution);
    mdf = GaussExc2MDF(exc, NA, n0, n, n1, focpos, lamem, mag, av, zpin);
    [modres, autotime] = FCS([exc.rho(:,1:end-1), exc.rho],[-exc.z(:,end:-1:2)+focpos exc.z-focpos],[mdf.volx(:,end:-1:2,:), mdf.volx],[mdf.voly(:,end:-1:2,:), mdf.voly]);

    u = squeeze(mdf.volx(:,:,1) + mdf.voly(:,:,1));
    close
    
	veff = 2*DetectionVolume(exc.rho,exc.z,mdf.volx,mdf.voly);
end

% physics
BoltzmannConstant = 1.380650324e-23*1e7;
AvogadroConstant = 6.0221419947e23;
Temperature = 273.15 + 36;
BeadDiam = 5e-2; % bead diam in mum
Eta = 1e-2; % solvent viscosity dyn/cm^2*s

% numerics
conc = 1e4; % number of particles per mum^3
N = 2^17; % max number of time channels
auto = zeros(N,1);
step = 0.01;
K = 1e2; % numer of markers

% generate particle distribution
x = -rhomax:step:rhomax;
y = -rhomax:step:rhomax;
z = -zmax:step:zmax;
%particles = poissrnd(conc*step^3,length(x), length(y), length(z));
cx = randi([2 length(x)-1],1,K);
cy = randi([2 length(y)-1],1,K);
cz = randi([2 length(z)-1],1,K);

% initialization
fluo = zeros(N,K);

% 2nd virial
al = 0.1;

for jt = 1:N
    for k=1:K
        xx = cx(k); xl = xx-1; xr = xx+1;
        yy = cy(k); yl = yy-1; yr = yy+1;
        zz = cz(k); zl = zz-1; zr = zz+1;
        if rand<(1+al*(particles(xx,yy,zz)-particles(xr,yy,zz)))/(2+al*(2*particles(xx,yy,zz)-particles(xr,yy,zz)-particles(xl,yy,zz)))
            cx(k) = cx(k)+1;
            if cx(k)==length(x)
                cx(k) = 2;
            end
        else
            cx(k) = cx(k)-1;
            if cx(k)==1
                cx(k) = length(x)-1;
            end
        end
        if rand<(1+al*(particles(xx,yy,zz)-particles(xx,yr,zz)))/(2+al*(2*particles(xx,yy,zz)-particles(xx,yr,zz)-particles(xx,yl,zz)))
            cy(k) = cy(k)+1;
            if cy(k)==length(y)
                cy(k) = 2;
            end
        else
            cy(k) = cy(k)-1;
            if cy(k)==1
                cy(k) = length(y)-1;
            end
        end
        if rand<(1+al*(particles(xx,yy,zz)-particles(xx,yy,zr)))/(2+al*(2*particles(xx,yy,zz)-particles(xx,yy,zr)-particles(xx,yy,zr)))
            cz(k) = cz(k)+1;
            if cz(k)==length(z)
                cz(k) = 2;
            end
        else
            cz(k) = cz(k)-1;
            if cz(k)==1
                cz(k) = length(z)-1;
            end
        end
    end
    tmp = particles.*(1 + al*particles);
    particles = (tmp([2:end  1],:,:) + tmp([end 1:end-1],:,:) + ...
        tmp(:,[2:end  1],:) + tmp(:,[end 1:end-1],:,:) + ...
        tmp(:,:,[2:end  1]) + tmp(:,:,[end 1:end-1]))/6 - tmp;
    
    % fluorescence
    for k=1:K
        rr = sqrt(x(cx(k)).^2 + y(cy(k)).^2);
        fluo(jt,k) = interp2(exc.z-focpos, exc.rho, u, abs(z(cz(k))), rr);
    end
    
    jt
    
end

flc = fft(sum(fluo,2));
flc = abs(ifft(flc.*conj(flc)));
auto = auto + flc(1:N);
semilogx(time(2:end/2), auto(2:end/2)); drawnow;

disp(k)

time = time(2:end/2);
auto = auto(2:end/2);
semilogx(time,auto/max(auto),autotime*5e-5/dif,(modres+K/8/rhomax^2/zmax*veff)/(1+K/8/rhomax^2/zmax*veff))
xlabel('time (s)'); ylabel('autocorrelation (a.u.)')
