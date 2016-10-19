% Atto DSOM

tau = 2;
Tpulse = 0.05;
Trep = 25;
kisc = 0;
extinction = 120e3/2;
wavelength = 640;
focus = 250;

x=1:200; 
y = PulsedExcitation(x,1/tau,Tpulse,Trep,kisc,extinction,wavelength,focus);

% saturation curve
plot(x,y*1e9/1e6)
xlabel('excitation power [\muW]'); ylabel('fluorescence emission rate [MHz]')
text(1/3,1/2,{'extinction @ 640 nm: 60,000 l/mol/cm','rep. rate: 40 MHz','focus diameter: 250 nm',...
    'quantum yield & detection efficiency: 1'},'units','normal')

% laser profile
rho = 0:1e3;
z = exp(-2*rho.^2/250^2);

% measured signal
yy = PulsedExcitation(z'*x,1/tau,Tpulse,Trep,kisc,extinction,wavelength,focus);
yd = rho*yy;

% master curves
ym = PulsedExcitation([1;0.5;0.1]*x,1/tau,Tpulse,Trep,kisc,extinction,wavelength,focus);
rhom = round(sqrt(-log([1;0.5;0.1])*focus^2/2));

% saturation vs. position
plot(x,ym./(diff(ym(:,1:2),[],2)*ones(1,size(ym,2))))
colorize
hold on
plot(x,yd/diff(yd(1:2)),'og')
hold off
xlabel('excitation power [\muW]')
ylabel('fluorescence emission rate [a.u.]')
legend({['\rho = ' int2str(rhom(1)) ' nm'],...
    ['\rho = ' int2str(rhom(2)) ' nm'],...
    ['\rho = ' int2str(rhom(3)) ' nm'],...
    'homoegeneous fluorophore distribution'},2)

% measurement regime
xx = 10:10:200;
tt = cumsum(1e3./xx/sum(1./xx));
bar(xx,tt)
axis([0 210 0 300])
xlabel('excitation power [\muW]'); ylabel('measurement time [\mus]')

% saturation vs. time
t = 0:1e3;
tmp = ones(size(t));
for j=2:length(tt) tmp = tmp + (t>tt(j-1)); end
x = xx(tmp);
ym = PulsedExcitation([1;0.5;0.1]*x,1/tau,Tpulse,Trep,kisc,extinction,wavelength,focus);
rhom = round(sqrt(-log([1;0.5;0.1])*focus^2/2));
yy = PulsedExcitation(z'*x,1/tau,Tpulse,Trep,kisc,extinction,wavelength,focus);
yd = rho*yy;

plot(t,ym./(ym(:,1)*ones(1,size(ym,2))))
colorize
hold on
plot(t,yd/yd(1),'og')
hold off
xlabel('time [\mus]')
ylabel('fluorescence emission rate [a.u.]')
legend({['\rho = ' int2str(rhom(1)) ' nm'],...
    ['\rho = ' int2str(rhom(2)) ' nm'],...
    ['\rho = ' int2str(rhom(3)) ' nm'],...
    'homoegeneous fluorophore distribution'},2)

