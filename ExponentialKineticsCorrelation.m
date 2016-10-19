state_vec = 1; % state vector of emitters
nmax=1;
tmax = 1e5;
kp = 0.01;
km = 0.001;

for k=1:tmax
    state_vec(k+1) = state_vec(k) - (state_vec(k)==1).*(rand(nmax,1)<=kp) + (state_vec(k)==0).*(rand(nmax,1)<=km);
end

