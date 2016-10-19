function p = DertiOptSat(w0,a0,pow,method)

global pd
pd =  1/5e-5/1e6;

natconst
extinction = 125e3/2;
wavelength = 640e-9;
tau = 2;
Tpulse = 0.05;
Trep = 50;
kisc = 0;
lin = 2*extinction*log(10)/AvogadroConstant*pow*1e-6/(PlanckConstant*SpeedOfLight/wavelength)/(pi*1e-18)/1e10;
lamex = 640/1.333;
lamem = 670/1.333;
av = 100/60*1e3;
close; 
rho = (0:0.025:5)*1e3;
if nargin<4 || isempty(method)
    z = (-5:0.5:5)*1e3;
    [z,rho]=meshgrid(z,rho);
    volx0 = exp(-2*rho.^2./LaserBeamFit(w0,lamex,z).^2)./LaserBeamFit(w0,lamex,z).^2.*D1CEF(a0,av,lamem,z);
    for j=1:length(pow)
        sat=pow(j);
        if sat>0
            volx = PulsedExcitation(lin(j)*volx0,1/tau,Tpulse,Trep,kisc);
            [modres, autotime] = FCS(1e-3*rho,1e-3*z,volx,volx);
            p(:,j) = simplex('Gauss2Fcs',[w0 w0 a0],[],[],[],[],av,[lamex lamem],[],autotime,modres(:,1),[],[],[pd;pd]);
        else
            p(:,j) = [w0; w0; a0];
        end
    end
elseif method==2
    volx0 = exp(-2*rho.^2./w0^2)./w0^2;
    for j=1:length(pow)
        sat=pow(j);
        if sat>0
            volx = PulsedExcitation(lin(j)*volx0,1/tau,Tpulse,Trep,kisc);
            tmp = simplex('Gauss',[0 w0],[0 0],[0 inf],[],[],rho,volx,[],[],[]);
            p(j) = 2*tmp(2);
        else
            p(j) = w0;
        end
    end
end
