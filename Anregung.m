%lambda = 635e-9;
lambda = 514e-9;
extinction = 1e5;
tau = 3;
rep = 1e3;
pulse = 1e3;
kisc = 1e2;

NatConst;
hnu = PlanckConstant*SpeedOfLight/lambda;

power = (0.001:0.001:1)*1e-3; % Anregungsleistung in W/cm^2
area = pi*(1e-4*0.2)^2; 
intensity = power/area;

sigma = extinction*log(10)*1e3/AvogadroConstant;

x = sigma*intensity/hnu/1e9;

plot(power*1e3,PulsedExcitation(x,1/tau,pulse,rep,kisc)/1e-9); figure(gcf)
xlabel('excitation power [mW]');
ylabel('fluorescence rate [1/s]');