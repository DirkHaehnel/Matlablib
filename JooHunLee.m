t = (0:512)*1e-10;
a = 1e3;
IRF = a*exp(-((t-30e-10)/2/5e-10).^2);
tau1 = 2*1e-9;
tau2 = 5*1e-9;
A1 = 1;
A2 = 1;
signal = A1*exp(-t/tau1)+A2*exp(-t/tau2); %TCSPC Count1
poi_signal = conv(signal,IRF); 
poi_signal = poi_signal(1:length(t));
poi_signal = poissrnd(poi_signal); %TCSPC Count2
close; plot(t,IRF,t,poi_signal)

[cx, taux, offset, csh, z] = DistFluofit(IRF, poi_signal, max(t), mean(diff(t)), [], 1);

