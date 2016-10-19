% program KaDe
clear all

kd = 2.3e-7; % K_D value
s0 = 10.^(-9:0.1:-3); % vector of initial SR concentrations
h0 = s0'; % vector of initial HTR concentrations
tsh = 0.9; % desired binding value

f = inline('((s0+x+kd)/2-sqrt((s0+x+kd).^2/4-s0*x)-tsh*s0).^2','x','s0','kd','tsh');

for j=1:length(s0)
    res(j) = fminbnd(f, 0, 1, optimset('TolX',1e-12), s0(j), kd, tsh);
end

loglog(s0,res,'r',s0,s0,'b',s0,sign(s0)*tsh*kd/(1-tsh),':g','linewidth',2); set(gca,'fontname','times','fontsize',16)
axis([s0(1) s0(end) kd max(res)])
xlabel('[SR_0]/M');
ylabel(['[HTR_0]/M for ' num2str(tsh*100) ' % binding'])
grid
