function [int, rho, phi, intc, ints] = QDSEP(kappa,psi,om,om0,maxrho,z0,NA,n0,n,n1,d0,d,d1,lambda,mag,focus,phi0);

% Quantum Dot Images

[tmp, tmp, tmp, rho, phi, fxx0, fxx2, fxz, byx0, byx2, byz] = SEPDipole([0 maxrho],z0,NA,n0,n,n1,d0,d,d1,lambda,mag,focus);

if nargin<17 | isempty(phi0)
    phi0 = 0;
end
phi = phi + phi0;

col = ones(size(phi));

f00 = fxx0.*conj(byx0);
f01 = -fxx0.*conj(byz);
f02 = -fxx0.*conj(byx2);
f10 = -fxz.*conj(byx0);
f11 = fxz.*conj(byz);
f12 = fxz.*conj(byx2);
f20 = -fxx2.*conj(byx0);
f21 = fxx2.*conj(byz);
f22 = fxx2.*conj(byx2);

intc(1,:) = 3*f00+3*f22+f11 - f00*cos(2*psi)*(1-cos(2*om)+kappa*(3+cos(2*om))*cos(2*om0)) + ...
    (f00-f11+f22)*(cos(2*om)+2*kappa*cos(2*om0)*sin(om)^2) + 4*kappa*f00*cos(om)*sin(2*psi)*sin(2*om0);

intc(2,:) = -2*(2*f01-f21+2*f10-f12)*sin(om)*(cos(psi)*cos(om)*(1-kappa*cos(2*om0))+kappa*sin(psi)*sin(2*om0));

intc(3,:) = -3*f20+f11-3*f02 + (f20+f02)*cos(2*psi)*(kappa*(3+cos(2*om))*cos(2*om0)+2*sin(om)^2) - ...
    (f20+f11+f02)*(cos(2*om)+2*kappa*cos(2*om0)*sin(om)^2) - 4*(f20+f02)*kappa*cos(om)*sin(2*psi)*sin(2*om0);

intc(4,:) = 2*(f21+f12)*sin(om)*(cos(psi)*cos(om)*(1-kappa*cos(2*om0))+kappa*sin(psi)*sin(2*om0));

intc(5,:) = -f22*(cos(2*psi)*(kappa*(3+cos(2*om))*cos(2*om0)+2*sin(om)^2)-4*kappa*cos(om)*sin(2*psi)*sin(2*om0));

ints(1,:) = 2*(f21+f12)*sin(om)*(cos(om)*(1-kappa*cos(2*om0))*sin(psi)-kappa*cos(psi)*sin(2*om0));

ints(2,:) = (f20+f02)*(sin(2*psi)*(kappa*(3+cos(2*om))*cos(2*om0)+2*sin(om)^2)+4*kappa*cos(2*psi)*cos(om)*sin(2*om0));

ints(3,:) = ints(1,:);

ints(4,:) = -f22*(sin(2*psi)*(kappa*(3+cos(2*om))*cos(2*om0)-2*sin(om)^2)+4*kappa*cos(2*psi)*cos(om)*sin(2*om0));

int = col*intc(1,:);
for j=1:4 int = int + cos(j*phi)*intc(j+1,:) + sin(j*phi)*ints(j,:); end
int = real(int);

