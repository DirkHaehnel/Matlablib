% Test formalism of FLORTRAN program as described in:
% Isenberg, I. et al., Biophys. J. 13 (1973), 1090 - 1115

N = 10000;

t_1 = 2;
t_2 = 4;
A1 = 0.50./t_1;
A2 = 0.75./t_2;

% Perfect cumulants:

C = [(t_1*A1 + t_2*A2) (t_1^2*A1 + t_2^2*A2) (t_1^3*A1 + t_2^3*A2) (t_1^4*A1 + t_2^4*A2)]./(A1 + A2);

fprintf('%5.3f %5.3f %5.2f %5.2f \n', C);
fprintf('%5.2f %5.2f %5.2f %5.3f \n', [t_1 A1/(A1+A2) t_2 A2/(A1+A2)]);

% Compute ideal TCSPC-histogram

Resolution = 0.008;
tau = (Resolution/2:Resolution:25);

IRF = exp(-(tau-1).^2./(0.08.^2));
dec = A1.*exp(-tau./t_1) + A2.*exp(-tau./t_2);
p = conv(IRF, dec);
p = p(1:numel(tau));

tcspc = p./sum(p).*N;

% Exponential depression

lambda = -log(0.5)/tau(end);
DEP    = exp(-lambda.*tau);

% Moments of IRF

H    = sum(          IRF.*DEP);
H(2) = sum((tau   ).*IRF.*DEP);
H(3) = sum((tau.^2).*IRF.*DEP);
H(4) = sum((tau.^3).*IRF.*DEP);

% Convert units to cnts/ns (see p. 1107)

H = H./Resolution;

% Moments of fluorescence signal

F    = sum(          tcspc.*DEP);
F(2) = sum((tau   ).*tcspc.*DEP);
F(3) = sum((tau.^2).*tcspc.*DEP);
F(4) = sum((tau.^3).*tcspc.*DEP);

% Estimate mean lefetime and amplitude (see p. 1107)

tau_m = F(2)/F(1) - H(2)/H(1);
amp_m = F(1)/H(1)/tau_m;

% Perform scaling of moments according to eq. A10 and A11

H(2) = H(2)/tau_m;
H(3) = H(3)/tau_m^2;
H(4) = H(4)/tau_m^3;

F(1) = F(1)/amp_m/tau_m;
F(2) = F(2)/amp_m/tau_m^2;
F(3) = F(3)/amp_m/tau_m^3;
F(4) = F(4)/amp_m/tau_m^4;

rep = 100;
while rep>0

% compute cumulants 

    G(1) =  F(1)                                           /H(1);
    G(2) = (F(2)   - G(1)*H(2)                            )/H(1);
    G(3) = (F(3)/2 - G(1)*H(3)/2 - G(2)*H(2)              )/H(1);
    G(4) = (F(4)/6 - G(1)*H(4)/6 - G(2)*H(3)/2 - G(3)*H(2))/H(1);

% compute lifetimes as roots of eq. 6

    k_1 = (G(2)*G(3) - G(1)*G(4))/(G(2)*G(2) - G(1)*G(3))/2;
    k_2 = (G(2)*G(4) - G(3)*G(3))/(G(2)*G(2) - G(1)*G(3));

    tau_1 = (k_1 - sqrt(k_1^2 + k_2));
    tau_2 = (k_1 + sqrt(k_1^2 + k_2));

% calculate according amplitudes

    amp_1 = (G(2) - G(1).*tau_2)./(tau_1.*(tau_1-tau_2));
    amp_2 = (G(2) - G(1).*tau_1)./(tau_2.*(tau_2-tau_1));

% cutoff correction

% beta is calculated according to eq. 10
    
    b(1) = amp_1 * sum(IRF.*exp(tau./tau_1))./Resolution;
    b(2) = amp_2 * sum(IRF.*exp(tau./tau_2))./Resolution;

% I is computed according to eq. 11

    I(1,1) =            tau_1*exp(-tau(end)/(tau_1));
    I(2,1) =            tau_2*exp(-tau(end)/(tau_2));
    I(1,2) = tau(end)  *tau_1*exp(-tau(end)/(tau_1))+   tau_1*I(1,1);
    I(2,2) = tau(end)  *tau_2*exp(-tau(end)/(tau_2))+   tau_2*I(2,1);
    I(1,3) = tau(end)^2*tau_1*exp(-tau(end)/(tau_1))+ 2*tau_1*I(1,2);
    I(2,3) = tau(end)^2*tau_2*exp(-tau(end)/(tau_2))+ 2*tau_2*I(2,2);
    I(1,4) = tau(end)^3*tau_1*exp(-tau(end)/(tau_1))+ 3*tau_1*I(1,3);
    I(2,4) = tau(end)^3*tau_2*exp(-tau(end)/(tau_2))+ 3*tau_2*I(2,3);
        
    delta = b*I;
    
    F = F + delta;    
    
    rep = rep -1;
end

fprintf('%5.2f %5.3f %5.2f %5.3f   ', [tau_1*tau_m amp_1/(amp_1+amp_2) tau_2*tau_m amp_2/(amp_1+amp_2)]);
fprintf('%5.3f %5.3f %5.2f %5.2f \n', G);

dec = amp_1*amp_m.*exp(-tau./tau_1/tau_m) + amp_2*amp_m.*exp(-tau./tau_2/tau_m);
zz = conv(IRF, dec);
zz = zz(1:numel(tau));
semilogy(tau,tcspc,tau,zz);
axis([0 tau(end) 1 10*ceil(N/numel(tau))]);


