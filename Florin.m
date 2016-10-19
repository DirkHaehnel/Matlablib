% Florin

% refractive indices
nw = 1.333; % water
ng = nw; % 1.51; % glass
nb = 1.5; %0.1937 + 3.2781i % bead
% wavelength
lambda = 0.633;
% NA
NA = 1;

% focpos = 0;

% coordinates of the bead interface
n = 300;
theta = (0.5:n)/n*pi;
r = ones(size(theta));

zv = r.*cos(theta);
rhov = r.*sin(theta);

maxnum = 1e3; % max number of plane waves for focus calculation
maxhar = 25; % max order of used spherical harmonics

kw = 2*pi/lambda*nw; kg = 2*pi/lambda*ng; kb = 2*pi/lambda*nb; 

% calculation of the excitation field at the bead surface
chimin = 0;
chimax = asin(NA/nw);
dchi = chimax/maxnum;
chi = (chimin:dchi:chimax)'; 

cc = cos(chi);
ss = sin(chi);
fac = ss.*sqrt(cc);

phase1 = kw*cc*(zv-focpos);
ez1 = exp(i*phase1);

barg = kw*ss*rhov;
j0 = besselj(0,barg); j1 = besselj(1,barg); j2 = besselj(2,barg);
fx0 = sum(((fac.*(cc+1))*ones(size(zv))).*(ez1.*j0)); % cos(0*phi)-component
fx2 = -sum(((fac.*(cc-1))*ones(size(zv))).*(ez1.*j2)); % cos(2*phi)-component
fz = -2*i*sum(((fac.*ss)*ones(size(zv))).*(ez1.*j1)); % cos(phi)-component

% er = (fx0+fx2).*sin(theta) + fz.*cos(theta); % cos(phi)
et = (fx0+fx2).*cos(theta) - fz.*sin(theta); % cos(phi)
ef = i*(fx0-fx2).*ones(size(theta)); % sin(phi)
% br = nw*er; % sin(phi)
bt = -i*nw*et; % sin(phi)
bf = -i*nw*ef; % cos(phi)

% calculation of the scattering matrix
% order of coefficients:
% M11b, N11b, M21b, N21b, ...M11w, N11w, M21w, N21w, ...,

r = r'; theta = theta';
At = []; Af = []; Bt = []; Bf = [];
for j=1:maxhar
    [mt, mf] = VSHM(j,1,kb*r,theta);
    [tmp, nt, nf] = VSHN(j,1,kb*r,theta);
    At = [At mt];
    Af = [Af mf];
    At = [At nt];
    Af = [Af nf];
    Bt = [Bt -i*nb*nt];
    Bf = [Bf -i*nb*nf];
    Bt = [Bt -i*nb*mt];
    Bf = [Bf -i*nb*mf];
end
for j=1:maxhar
    [mt, mf] = VSHM(j,1,kw*r,theta,'h');
    [tmp, nt, nf] = VSHN(j,1,kw*r,theta,'h');
    At = [At -mt];
    Af = [Af -mf];
    At = [At -nt];
    Af = [Af -nf];
    Bt = [Bt i*nw*nt];
    Bf = [Bf i*nw*nf];
    Bt = [Bt i*nw*mt];
    Bf = [Bf i*nw*mf];
end

c = [At; Af; Bt; Bf]\[et.'; ef.'; bt.'; bf.'];

At = []; Af = []; Bt = []; Bf = [];
for j=1:maxhar
    [mt, mf] = VSHM(j,1,kw*r,theta);
    [tmp, nt, nf] = VSHN(j,1,kw*r,theta);
    At = [At mt];
    Af = [Af mf];
    At = [At nt];
    Af = [Af nf];
    Bt = [Bt -i*nw*nt];
    Bf = [Bf -i*nw*nf];
    Bt = [Bt -i*nw*mt];
    Bf = [Bf -i*nw*mf];
end

c0 = [At; Af; Bt; Bf]\[et.'; ef.'; bt.'; bf.'];

[aat,aaf] = VSH2ADR(1:maxhar,ones(1,maxhar),zeros(1,maxhar),0);
[bbt,bbf] = VSH2ADR(1:maxhar,ones(1,maxhar),ones(1,maxhar),0);
st = [st aat*c((1:2:2*maxhar)+2*maxhar) + bbt*c((2:2:2*maxhar)+2*maxhar)];
sf = [sf aaf*c((1:2:2*maxhar)+2*maxhar) + bbf*c((2:2:2*maxhar)+2*maxhar)];
st0 = [st0 aat*c(1:2:2*maxhar) + bbt*c(2:2:2*maxhar)];
sf0 = [sf0 aaf*c(1:2:2*maxhar) + bbf*c(2:2:2*maxhar)];

% r1 = (0.5:50)/50;
% r2 = max(r1)+r1;
% t1 = zeros(length(theta),length(r1)); t2 = t1; f1 = t1; f2 = t1;
% for j=1:maxhar 
%     [aat, aaf] = VSHM(j,1,kb*r1,theta);
%     [bbr, bbt, bbf] = VSHN(j,1,kb*r1,theta);
%     t1 = t1 + c(2*j-1)*(cos(theta)*ones(size(r1))).*aat + ...
%         c(2*j)*((cos(theta)*ones(size(r1))).*bbt+(sin(theta)*ones(size(r1))).*bbr); 
%     f1 = f1 + c(2*j-1)*aaf + c(2*j)*bbf;     
%     [aat, aaf] = VSHM(j,1,kw*r2,theta);
%     [bbr, bbt, bbf] = VSHN(j,1,kw*r2,theta);
%     t2 = t2 + c(2*(maxhar+j)-1)*(cos(theta)*ones(size(r2))).*aat + ...
%         c(2*(maxhar+j))*((cos(theta)*ones(size(r2))).*bbt+(sin(theta)*ones(size(r2))).*bbr); 
%     f2 = f2 + c(2*(maxhar+j)-1)*aaf + c(2*(maxhar+j))*bbf;  
% end
% f1 = i*f1;
% f2 = i*f2;
% 
% rr = [r1 r2];
% 
% % barg = kw*rr'*ss';
% % j0 = besselj(0,barg); j1 = besselj(1,barg); j2 = besselj(2,barg);
% % ffx0 = j0*(((fac.*(cc+1))*ones(size(zv))).*ez1); % cos(0*phi)-component
% % ffx2 = j2*(-((fac.*(cc-1))*ones(size(zv))).*ez1); % cos(2*phi)-component
% % ffz = 2*i*j1*(((fac.*ss)*ones(size(zv))).*ez1); % cos(phi)-component
% 
% close all
% pcolor(sin(theta)*[-fliplr(rr) rr],cos(theta)*[fliplr(rr) rr],[fliplr(real([t1 t2])) real([t1 t2])]); shading interp; axis image
% drawnow
% % figure; pcolor([-fliplr(rr) rr],zv,[flipud(real(ffx0+ffx2)); real(ffx0+ffx2)]'); shading interp; axis image
