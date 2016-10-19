function [kappa1,kappa2,bin] = DICTriplet(exc,mdf,pulse,sat,triplet)

maxm = exc.maxm;
rho = exc.rho;
z = exc.z;
fxc1 = exc.fxc1;
fxs1 = exc.fxs1;
fyc1 = exc.fyc1;
fys1 = exc.fys1;
fzc1 = exc.fzc1;
fzs1 = exc.fzs1;
fxc2 = exc.fxc2;
fxs2 = exc.fxs2;
fyc2 = exc.fyc2;
fys2 = exc.fys2;
fzc2 = exc.fzc2;
fzs2 = exc.fzs2;

volx1 = mdf.volx1;
voly1 = mdf.voly1;
volx2 = mdf.volx2;
voly2 = mdf.voly2;

if nargin<3 | isempty(pulse)
    pulse = [1 1];
end
if length(pulse)==1
    pulse = [pulse pulse];
end

if nargin<4 | isempty(sat)
    sat = 0;
end

if nargin<5 | isempty(triplet)
    triplet = 0;
end

tmp = 0; % following procedure is only approximative
for j=0:maxm
    psi = j/(2*maxm+1)*2*pi;
    tmp = max([tmp, max(max(FocusImage2D(exc.rho,exc.z,cat(4,cat(3,exc.fxc1,exc.fxs1),cat(3,exc.fyc1,exc.fys1),cat(3,exc.fzc1,exc.fzs1)),psi,'calc')))]);
end
sat1 = sat/tmp;
tmp = 0;
for j=0:maxm
    psi = j/(2*maxm+1)*2*pi;
    tmp = max([tmp, max(max(FocusImage2D(exc.rho,exc.z,cat(4,cat(3,exc.fxc2,exc.fxs2),cat(3,exc.fyc2,exc.fys2),cat(3,exc.fzc2,exc.fzs2)),psi,'calc')))]);
end
sat2 = sat/tmp;
bin = 2*(0.5:100)/100*PulsedExcitation(sat,1,pulse(1),pulse(2));
for j=0:2*maxm
    psi = j/(2*maxm+1)*2*pi;
    vx1 = FocusImage2D(exc.rho,exc.z,mdf.volx1+mdf.voly1,psi,'calc');
    vx2 = FocusImage2D(exc.rho,exc.z,mdf.volx2+mdf.voly2,psi,'calc');
    int = PulsedExcitation(sat1*FocusImage2D(exc.rho,exc.z,cat(4,cat(3,exc.fxc1,exc.fxs1),cat(3,exc.fyc1,exc.fys1),cat(3,exc.fzc1,exc.fzs1)),psi,'calc'),1,pulse(1),pulse(2))+...
        PulsedExcitation(sat2*FocusImage2D(exc.rho,exc.z,cat(4,cat(3,exc.fxc2,exc.fxs2),cat(3,exc.fyc2,exc.fys2),cat(3,exc.fzc2,exc.fzs2)),psi,'calc'),1,pulse(1),pulse(2));
  
    kappa1 = mhist(int(:),bin,int(:).*vx1(:).*exc.rho(:));
    kappa2 = mhist(int(:),bin,int(:).*vx2(:).*exc.rho(:));    

end

