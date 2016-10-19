function [y z] = ExcitonExcitation(x,k10,extinction,Tpulse,Trep,wavelength)

AvogadroConstant = 6.0221419947e23; % NIST
PlanckConstant = 6.6260687652e-34; % NIST
SpeedOfLight = 299792458; % exact
dt = 1e-3; % integration time unit

if nargin <6 || isempty(wavelength)
    wavelength = 470; % nm
end
% x in muW/mum^2
sig = extinction*log(10)/AvogadroConstant/(PlanckConstant*SpeedOfLight/wavelength/1e-9)*1e-2;  % photons/ns!!!

x = x*(Trep+Tpulse)/Tpulse;
y = zeros(length(k10)+1,size(x,1),size(x,2),size(x,3));
z = zeros(length(k10),size(x,1),size(x,2),size(x,3));

M1 = zeros(length(k10)+1); M2 = M1;
M1(1,1) = -sig(1);
M2(1,2:end) = k10(:)';
for j=2:length(k10)+1
    M1(j,j-1) = sig(j-1);
    if j<=length(k10)
        M1(j,j) = -sig(j);
    end
    M2(j,j) = -k10(j-1);
end
for jx=1:size(x,1)
    for jy=1:size(x,2)
        for jz=1:size(x,3)
            M = expm((Trep-Tpulse)*M2)*expm(Tpulse*(x(jx,jy,jz)*M1+M2));
            old = [1 0 0 0]';
            tst=M*old;
            while sum(abs(tst-old).^2)>1e-10
                old=tst;
                tst=M*old;
            end
            y(:,jx,jy,jz) = tst;

            M = expm(dt*(x(jx,jy,jz)*M1+M2));
            for k=1:ceil(Tpulse/dt)
                z(:,jx,jy,jz) = z(:,jx,jy,jz) + diag(k10)*dt*tst(2:end);
                tst = M*tst;
            end
            z(:,jx,jy,jz) = 1e9/(Trep+Tpulse)*(z(:,jx,jy,jz) - (expm(M2(2:end,2:end)*(Trep-Tpulse))-eye(length(k10)))*tst(2:end));
        end
    end
end
y = squeeze(y);
z = squeeze(z);

return

[yc zc] = ExcitionExcitation(x,[1/20 1/0.7 1/0.5],2e5*[1 1 1],0.5,0.5);
[yp zp] = ExcitionExcitation(x,[1/20 1/0.7 1/0.5],2e5*[1 1 1],0.05,12.5);

k = 2;
loglog(x,zp(k,:),x,zc(k,:),x,zp(k,1)/x(1).^k*x.^k,'r:',x,zc(k,1)/x(1).^k*x.^k,'b:')