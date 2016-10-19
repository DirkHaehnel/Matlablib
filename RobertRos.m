% program Robert Ros

clear all 
close all

lamex = 488;
lamem = 585;
nair = 1.;
nglass = 1.52;
nsiex = 3.9 + i*0.017;
nsiem = 3.9 + i*0.017;
load metals
nauex = gold(wavelength==lamex);
nauem = gold(wavelength==lamem);
dau = 20; % thickness of silver layer in nm

NA = 1.45;

z = 0;
d = 20:2:1000;
theta = (0.5:5e2)/5e2*asin(NA/nglass);
st = sin(theta)*ones(1,length(z));
clear tirfp tirfv vk pk vt pt res

incidence = 42:2:60;
for k=1:length(incidence)
    % excitation
    w0 = nglass*cos(incidence(k)/180*pi);
    w1 = sqrt(nair^2-nglass^2+w0^2);
    [bla,bla,tp] = fresnel(w0,nglass,nair);
    rp = fresnel(w1,nair,nglass);
    rp1 = fresnel(w1,[nair,nauex,nsiex],dau/lamex*2*pi);
    for j=1:length(d)
        tirfv(j,k) = abs(1-(w1./nair).^2)*abs(tp*(1-rp1*exp(i*2*pi*d(j)/lamex*w1))/(1-rp*rp1*exp(i*2*pi*d(j)/lamex*w1))).^2;
        tirfp(j,k) = abs(w1./nair).^2*abs(tp*(1-rp1*exp(i*2*pi*d(j)/lamex*w1))/(1-rp*rp1*exp(i*2*pi*d(j)/lamex*w1))).^2;        
    end
end
for j=1:length(d)
    [v,pc,ps] = DipoleL(theta,0,nglass,nair,[nauem nsiem],[],d(j)/lamem*2*pi,dau/lamem*2*pi);
    [bla,bla,bla,bla,qvd,qvu,qpd,qpu] = LifetimeL(0,nglass,nair,[nauem nsiem],[],d(j)/lamem*2*pi,dau/lamem*2*pi);
    vk(j) = st*abs(v).^2;
    pk(j) = st*(abs(pc).^2+abs(ps).^2)/2;
    vt(j) = qvu+qvd;
    pt(j) = qpu+qpd;
end
res = diff(theta(1:2))*real((vk-pk)./(vt-pt)+(vt.*pk-pt.*vk).*atan(sqrt((vt-pt)./pt))./sqrt(pt)./(vt-pt).^(3/2))';

fs = 20;
semilogy(d,((pk./pt)'*ones(1,size(tirfp,2))).*tirfp)
xlabel('air gap width (nm)','fontname','arial','fontsize',fs); 
ylabel('observable fluorescence (a.u.)','fontname','arial','fontsize',fs)
title('parallel dipole orientation','fontname','arial','fontsize',fs)
set(gca,'fontname','arial','fontsize',fs,'TickDir','out')
colorize
text(0.0185,0.095,{'incidence angle from 42° (blue) through 60° (red)','','in steps of 2°, \lambda_{ex} = 488 nm, \lambda_{em} = 585 nm'},'units','normalized','fontname','arial','fontsize',18)

return
semilogy(d,((vk./vt)'*ones(1,size(tirfv,2))).*tirfv)
xlabel('air gap width (nm)','fontname','arial','fontsize',fs); 
ylabel('observable fluorescence (a.u.)','fontname','arial','fontsize',fs)
title('vertical dipole orientation','fontname','arial','fontsize',fs)
set(gca,'fontname','arial','fontsize',fs,'TickDir','out')
colorize
text(0.0185,0.095,{'incidence angle from 42° (blue) through 60° (red)','','in steps of 2°, \lambda_{ex} = 488 nm, \lambda_{em} = 585 nm'},'units','normalized','fontname','arial','fontsize',18)




