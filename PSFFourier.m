function [ff, volx, rho, z] = PSFFourier(w0,a0,pow, rhofield, zfield)

natconst
extinction = 125e3/2;
wavelength = 640e-9;
tau = 2;
Tpulse = 0.05;
Trep = 50;
kisc = 0;
lin = 2*pow*extinction*log(10)/AvogadroConstant*pow*1e-6/(PlanckConstant*SpeedOfLight/wavelength)/(pi*1e-18)/1e10;
lamex = 640/1.333;
lamem = 670/1.333;
av = 100/60*1e3;
if nargin<4 || isempty(rhofield)
    rho = (0:0.025:5)*1e3;
else
    rho = (0:0.025:rhofield)*1e3;
end
if nargin<5 || isempty(zfield)
    z = (-10:0.025:10)*1e3;
else
    z = (-zfield:0.025:zfield)*1e3;
end
[z,rho]=meshgrid(z,rho);
volx = exp(-2*rho.^2./LaserBeamFit(w0,lamex,z).^2)./LaserBeamFit(w0,lamex,z).^2.*D1CEF(a0,av,lamem,z);
volx = PulsedExcitation(lin*volx,1/tau,Tpulse,Trep,kisc);
volx = [flipud(volx);volx];
rho = [-flipud(rho); rho]; z = [z;z];
ff = fftshift(abs(fft2(volx)).^2);


