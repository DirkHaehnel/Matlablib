N = 1e4; % number of photons

A1  = 1; % Amplitude of lifetime 1
t_1 = 3; % lifetime 1

C = [1 t_1  t_1^2 t_1^3 t_1^4]; % expectation values for 1st to 4th moments


Resolution = 0.016;
tau = (Resolution/2:Resolution:50);

IRF = exp(-(tau-1).^2./(0.08.^2));
dec = A1.*exp(-tau./t_1)./(t_1);    % decay function
p = conv(IRF, dec);                 % ideal tcspc histogram   
p = p(1:numel(tau));

tcspc = p./sum(p).*N;               % fill with according number of photons

plot(tau, tcspc)

H    = sum(         IRF);
H(2) = sum((tau   ).*IRF);
H(3) = sum((tau.^2).*IRF);
H(4) = sum((tau.^3).*IRF);
H(5) = sum((tau.^4).*IRF);

F    = sum(         tcspc);
F(2) = sum((tau   ).*tcspc);
F(3) = sum((tau.^2).*tcspc);
F(4) = sum((tau.^3).*tcspc);
F(5) = sum((tau.^4).*tcspc);

H = H/H(1);
F = F/F(1);

tav = 1;
for j=1:4
    tav(j+1) = 1/factorial(j)*F(j+1);
    for s=0:j-1
        tav(j+1) = tav(j+1) - 1/factorial(j-s)*H(j-s+1)*tav(s+1);
    end
end

display([C; tav]);
