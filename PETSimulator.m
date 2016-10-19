tmean = 100; % average time between photon detections in the bright state
kp = 0.01; % on-rate 
km = 0.01; % off-rate
N = 1e7; % ~ number of photons
N = round(N*(kp+km)/kp); 
photons = exprnd(tmean,[1 N]);

state = ones(1,N+1);
for j=1:N
    if state(j)==1
        state(j+1) = rand<kp/(kp+km)+km/(kp+km)*exp(-(kp+km)*photons(j));
    else
        state(j+1) = rand<kp/(kp+km)*(1-exp(-(kp+km)*photons(j)));
    end
end
t = cumsum([0 photons]);
t = round(t(state==1));    

return

[auto, autotime] = tttr2xfcs(t,ones(size(t)),10,10);
auto = auto/(max(t)-min(t));
semilogx(autotime,auto)

close; p=Simplex('ExpFun',156,[],[],[],[],autotime, auto, 3)

[~,c] = ExpFun(p,autotime, auto, 3)

[c(1) (1/tmean*kp/(kp+km))^2]
[c(2) 1/tmean^2*kp*km/(kp+km)^2]