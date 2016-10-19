function [y, y1, y2] = PulsedExcitation(x,k10,Tpulse,Trep,kisc,extinction,wavelength,w0)

% x - excitation chance per pulse or excitation power in muW
% k10 - inverse lifetime in 1/ns
% Tpulse - pulse duration in ns
% Trep - excitation period in ns
% kisc - ratio of intersystem crossing rate to phosphorescence rate
% extinction - in L/mol/cm
% wavelength - in nm
% w0 - focus diameter in nm
%
% y in 1/ns (or per pulse)

AvogadroConstant = 6.0221419947e23; % NIST
PlanckConstant = 6.6260687652e-34; % NIST
SpeedOfLight = 299792458; % exact

if nargin<5 || isempty(kisc)
    kisc = 0;
end
if nargin>5
    if isempty(extinction)
        extinction = 100e3; % l/mol/cm
    end
    if nargin<7 || isempty(wavelength)
        wavelength = 640; % nm
    end
    if nargin<8 || isempty(w0)
        w0 = 250; % nm
    end
    % x in muW
    x = 2*extinction*log(10)/AvogadroConstant*x*1e-6/(PlanckConstant*SpeedOfLight/wavelength/1e-9)/(pi*w0^2*1e-18)/1e10;  % photons/ns!!!
end

k01 = x*Trep/Tpulse;
% chance to be in the excited state S1: 
kappa = k01./(k01+k10)*Tpulse/Trep+k01.^2./(k01+k10).^2/Trep/k10.*(1-exp(-(k01+k10)*Tpulse)).*(1-exp(-k10*(Trep-Tpulse)))./(1-exp(-k10*Trep-k01*Tpulse));
y = k10*kappa./(1+kisc*kappa);

if nargout>1
kappa1 = (2.*exp(k10.*Trep).*k01.*k10 - 2.*exp(2.*(k01 + k10).*Tpulse + k10.*Trep).*k01.*k10 - 2.*exp(k01.*Tpulse + k10.*(Tpulse + Trep)).*(k01 + k10).*...
    (k01.^2 + k10^2).*Tpulse + exp(k01.*Tpulse + 2.*k10.*Trep).*k01.*(-2.*k10 + k01.*(k01 + k10).*Tpulse) + ...
    exp((k01 + 2.*k10).*Tpulse).*k01.*(2.*k10 + k01.*(k01 + k10).*Tpulse) + exp(k10.*Tpulse).*k10.*(-2.*k01 + k10.*(k01 + k10).*Tpulse) + ...
    exp(2.*k01.*Tpulse + k10.*Tpulse + 2.*k10.*Trep).*k10.*(2.*k01 + k10.*(k01 + k10).*Tpulse))./...
    (exp(k10.*Tpulse).*(-1 + exp(k01.*Tpulse + k10.*Trep)).^2.*k10.*(k01 + k10).^3.*Trep);
 
    y1 = k10./(1+kisc*kappa).^2.*kappa1.*Trep/Tpulse;
    
    if nargout>2
        kappa2 = -(((-6.*(-1 + exp((k01 + k10).*Tpulse)).*(-1 + exp(k10.*(-Tpulse + Trep))).*k01.^2)./...
            (-1 + exp(k01.*Tpulse + k10.*Trep)) + 2.*k10.^2.*(k01 + k10).*Tpulse - ((-1 + exp(k10.*(-Tpulse + Trep))).*...
            (k01 + k10).^2.*(-2 + 2.*exp(3.*k01.*Tpulse + k10.*Tpulse + 2.*k10.*Trep) + ...
            exp(2.*k01.*Tpulse + k10.*(Tpulse + Trep)).*(-4 + k01.*Tpulse.*(-4 + k01.*Tpulse)) - ...
            exp(2.*k01.*Tpulse + 2.*k10.*Trep).*(2 + k01.*Tpulse.*(-4 + k01.*Tpulse)) + ...
            exp(k01.*Tpulse + k10.*Trep).*(4 - k01.*Tpulse.*(4 + k01.*Tpulse)) + exp((k01 + k10).*Tpulse).*...
            (2 + k01.*Tpulse.*(4 + k01.*Tpulse))))./(-1 + exp(k01.*Tpulse + k10.*Trep)).^3 + ...
            (4.*(-1 + exp(k10.*(-Tpulse + Trep))).*k01.*(k01 + k10).*(2 + exp(k01.*Tpulse).*(-(exp(k10.*Tpulse).*...
            (2 + k01.*Tpulse)) + exp(k10.*Trep).*(-2 + 2.*exp((k01 + k10).*Tpulse) + k01.*Tpulse))))./...
            (-1 + exp(k01.*Tpulse + k10.*Trep)).^2)./(k10.*(k01 + k10).^4.*Trep));
        
        y2 = (-k10*kisc./(1+kisc*kappa).^3.*kappa1.^2 + k10./(1+kisc*kappa).^2.*kappa2)*(Trep/Tpulse)^2;
    end
end

return

if 0
    tf = 3;
    kisc = 1;
    x=0:0.01:1;
    tau=[0.05 1 5 12.5];
    for j=1:4 p(j)=simplex('OptSatMin',1,0,[],[],[],[-fliplr(x) x],[-PulsedExcitation(fliplr(x),1/tf,tau(j),25,kisc) PulsedExcitation(x,1/tf,tau(j),25,kisc)]); [bla,z(:,j)]=OptSatMin(p(j),x,PulsedExcitation(x,1/tf,tau(j),25,kisc)); end
    h = plot(x,PulsedExcitation(x,1/tf,0.05,25,kisc),x,PulsedExcitation(x,1/tf,1,25,kisc),x,PulsedExcitation(x,1/tf,5,25,kisc),x,PulsedExcitation(x,1/tf,12.5,25,kisc),x,PulsedExcitation(x,1/tf,25,25,kisc),x,z,'o','MarkerSize',4,'Linewidth',2);
    tmp = jet;
    for j=1:length(h) if j<=5 set(h(j),'color',tmp(12*j,:)); else set(h(j),'color',tmp(12*(j-5),:)); end; end
    times16
    axis tight
    xlabel('average excitation rate [GHz]');
    ylabel('fluorescence rate [GHz]');
    for j=1:length(tau)
        s{j} = [mnum2str(tau(j),3) ' ns']; 
    end
    s{j+1} = 'cw';
    legend(s,2)
end

if 0 % excitation rate
    NatConst
    extinction = 100e3;
    power = (1:3e2)*1e-6;
    wavelength = 400e-9;
    focus = 220e-9;
    sigma = 1e3*extinction*log(10)/AvogadroConstant;
    kex = sigma*power/(PlanckConstant*SpeedOfLight/wavelength)/(pi*focus^2)/1e4;
    tau = 2;
    plot(1e6*power,kex,1e6*power,1/tau/1e-9*(power>0))
    xlabel('power [\muW]'); ylabel('excitation rate [1/s]')

    % STED Lifetime
    plot(1e6*power,1./(kex+1/tau/1e-9))
    xlabel('power [\muW]'); ylabel('lifetime [ns]');
    
    tau = 4;
    Tpulse = 2;
    Trep = 25;
    kisc = 0;
    plot([0 1e6*power],[0 1e3*PulsedExcitation(kex/1e9,1/tau,Tpulse,Trep,kisc)],[0 1e6*power],[0 1e3*PulsedExcitation(kex/1e9,1/tau,Trep,Trep,kisc)],':')
    xlabel('laser excitation power [\muW]');
    ylabel('maximum excitation rate [MHz]')
    ax = axis;
    text(ax(1)+0.05*diff(ax(1:2)),ax(3)+0.8*diff(ax(3:4)),{['lifetime = ' mnum2str(tau,1,1) ' ns'],['1/e^2 focus radius = ' int2str(1e9*focus) ' nm'],...
        ['extinction = ' num2str(extinction) ' cm^{-1}\cdotM^{-1}'],['wavelength = ' int2str(1e9*wavelength) ' nm'],...
        ['repetition rate = ' mnum2str(1e3/Trep,2,1) ' MHz'],['pulse duration = ' int2str(1e3*Tpulse) ' ps']})
end

if 0
    natconst
    extinction = 100e3;
    power = (1:3e2)*1e-6;
    wavelength = 635e-9;
    focus = 200e-9;
    tau = 2;
    Tpulse = 0.05;
    Trep = 25;
    kisc = 0;
    
    lin = PulsedExcitation(extinction*log(10)/AvogadroConstant*power/(PlanckConstant*SpeedOfLight/wavelength)/(pi*focus^2)/1e10,1/0.1,Trep,Trep,kisc);
    pul = PulsedExcitation(extinction*log(10)/AvogadroConstant*power/(PlanckConstant*SpeedOfLight/wavelength)/(pi*focus^2)/1e10,1/tau,Tpulse,Trep,kisc);
    
    plot(power,PulsedExcitation(extinction*log(10)/AvogadroConstant*power/(PlanckConstant*SpeedOfLight/wavelength)/(pi*focus^2)/1e10,1/tau,Tpulse,Trep,kisc))
end

