% Programm zur Simulation eines FCS-Experiments

close all

fout = fopen('temp','w');
a = 0.25; % transversal PSF axis in mum
b = 2; % longitudinal PSF axis in mum
dt = 2e-7; % time interval in seconds
dif = [4e-6 2e-6]; % diffusion constants in cm^2/s unit
dif = dif*1e8;
rhomax = 8*a;
zmax = 8*b;
conc = [5e-10 1e-9]; % concentrations in M
AvogadroConstant = 6.0221419947e23;
M = round(8*rhomax^2*zmax*AvogadroConstant*1e-15*conc); % number of molecules
MM = sum(M);
xi = ones(1,MM);
M = cumsum([0 M]);
for kk=1:length(M)-1
    xi(M(kk)+1:M(kk+1)) = sqrt(2*dif(kk)*dt);
end
xi = diag(xi);
N = round(1e7/MM); % max number of time channels
L = (1:N)';
background = 0;
fluomax = dt*[1e6 5e5]; % max fluorescence rates for molecule in focus center
x = 2*(rand(1,MM)-0.5)*rhomax;
y = 2*(rand(1,MM)-0.5)*rhomax;
z = 2*(rand(1,MM)-0.5)*zmax;
offset = 0;
maxk = round(60/dt/N); % 60 seconds measurement time
for k=1:maxk 
   % random walk
   x = cumsum([x(end,:); randn(N,MM)*xi]);
   y = cumsum([y(end,:); randn(N,MM)*xi]);    
   z = cumsum([z(end,:); randn(N,MM)*xi]);
   x(1,:) = []; y(1,:) = []; z(1,:) = [];
   x = mod(x+rhomax,2*rhomax)-rhomax;
   y = mod(y+rhomax,2*rhomax)-rhomax;   
   z = mod(z+zmax,2*zmax)-zmax;   
   
   % fluorescence
   chan = []; tmp = [];
   for kk=1:length(M)-1
       ind = M(kk)+1:M(kk+1);
       fluo = sum(fluomax(kk)*exp(-2*(x(:,ind).^2+y(:,ind).^2)/a^2 - 2*z(:,ind).^2/b^2),2) + dt*background;
       fluo = rand(size(fluo))<1-exp(-fluo);
       chan = [chan; kk*fluo(fluo)];
       tmp = [tmp; offset + L(fluo)];
   end
   [tmp, ind] = sort(tmp);
   chan = chan(ind);
   fwrite(fout,[tmp chan]','uint32'); 
   offset = offset + N;
   disp(k/maxk)
end
fclose(fout);

fin = fopen('temp','r');
k = 1;
clear auto
while ~feof(fin)
    tttr = fread(fin,2e6,'uint32');
    tmp = tttr(2:2:end);
    tttr = tttr(1:2:end);
    chan = zeros(length(tmp),max(tmp(:)));
    for j=1:max(tmp(:))
        chan(:,j) = tmp==j;
    end
    if k==1
        [auto, autotime] = tttr2xfcs(tttr,chan,20,10);
        trace = tttr2bin(tttr, 40e-6/dt, tmp);
    else
        auto = auto + tttr2xfcs(tttr,chan,20,10);
        trace = [trace; tttr2bin(tttr, 40e-6/dt, tmp)];
    end
    k = k+1;
end
fclose(fin)


%[auto,autotime]=tttr2xfcs(tttr,ones(size(tttr)),20,10);
%p=Simplex('Rigler',[1e-4 1e-3], [], [], [], [], dt*autotime, auto, [], 1);