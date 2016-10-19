% comments
% don't trust the w1 w2 etc. values - they were fitted assuming that the
% Gauss approximation is valid!

load D:\MATLAB\2fFCSRadSatRef_full2.mat
satv = [0 0.05*exp(log(1/0.05)*(0:9)/9)];

% relation between laser beam diameter and MDF waist
rr = radv(1):radv(end);
p1 = Simplex('ExpFun',100,0,[],[],[],radv,squeeze(w1(:,1,1,1)),1);
[err,c1] = ExpFun(p1,radv,squeeze(w1(:,1,1,1)),1);
plot(rr/1e3,(c1(1)+c1(2)*exp(-(rr-min(rr))/p1))*1e3)
% p2 = Simplex('ExpFun',100,0,[],[],[],radv,squeeze(w2(:,1,1,1)),1);
% [err,c2] = ExpFun(p2,radv,squeeze(w2(:,1,1,1)),1);
% plot(rr/1e3,(c1(1)+c1(2)*exp(-(rr-min(rr))/p1))*1e3,rr/1e3,(c2(1)+c2(2)*exp(-(rr-min(rr))/p2))*1e3)
xlabel('laser beam radius [mm]'); 
ylabel('MDF waist [\mum]');

% dependence of fitted D on opt. saturation
plot(satv,squeeze(dc(:,2,1,:)))
xlabel('\itI_{max}\rm/\itI_{sat}')
ylabel('\itD_{fit}\rm/\itD')
for j=1:length(radv)
    text(0.8,interp1(satv,squeeze(dc(j,2,1,:)),0.8)+0.005,[mnum2str(radv(j)/1e3,1,2) ' mm'])
end
colorize

% dependence of fitted D on opt. sat. at 10 mum coverslide 
dc(1,end,1,3) = dc(1,end,2,3);
plot(satv,squeeze(dc(:,end,1,:)))
xlabel('\itI_{max}\rm/\itI_{sat}')
ylabel('\itD_{fit}\rm/\itD')
for j=1:length(radv)
    text(0.8,interp1(satv,squeeze(dc(j,end,1,:)),0.8)+0.005,[mnum2str(radv(j)/1e3,1,2) ' mm'])
end
colorize

% shift of autocorrelation upon opt. saturation
jrad = 3;
semilogx(autotime,squeeze(fullres(:,1,jrad,1,1,:))./(ones(size(fullres,1),1)*squeeze(fullres(1,1,jrad,1,1,:))'))
xlabel('time [a.u.]')
ylabel('autocorrelation')
colorize

% dependence of fitted D on opt. coverslide
dc(end,1,1,1) = dc(end,1,2,1);
plot(covv,squeeze(dc(:,:,1,1)))
xlabel('\delta [\mum]')
ylabel('\itD_{fit}\rm/\itD')
for j=1:length(radv)
    text(0.25,interp1(covv,squeeze(dc(j,:,1,1)),0.25)+0.003,[mnum2str(radv(j)/1e3,1,2) ' mm'])
end
colorize

% dependence of fitted D on opt. coverslide at max opt. sat.
dc(1,3,1,end) = dc(1,3,2,end);
dc(3,4,1,end) = dc(3,4,2,end);
plot(covv,squeeze(dc(:,:,1,end)))
xlabel('\delta [\mum]')
ylabel('\itD_{fit}\rm/\itD')
for j=1:length(radv)
    text(0.25,interp1(covv,squeeze(dc(j,:,1,end)),0.25)+0.005,[mnum2str(radv(j)/1e3,1,2) ' mm'])
end
colorize

% shift of autocorrelation upon coverslide
jrad = 1;
semilogx(autotime,squeeze(fullres(:,1,jrad,:,1,1))./(ones(size(fullres,1),1)*squeeze(fullres(1,1,jrad,:,1,1))'))
xlabel('time [a.u.]')
ylabel('autocorrelation')
colorize

% check of Gaussian assumption
lamex = 0.64; lamem = 0.67; n = 1.333; mag = 60; av = 100;
jrad = 5;
subplot(2,3,1:2)
plot(z-15,sqrt((wx1(:,jrad,1,1,1).^2+wy1(:,jrad,1,1,1).^2)/2),'o',z-15,LaserBeamFit(w1(jrad,1,2,1),lamex/n,z-15),'--',z-15,LaserBeamFit(w0(jrad,1,2,1),lamex/n,z-15));
xlabel('\itz\rm [\mum]')
ylabel('\itw\rm(\itz\rm) [\mum]')
subplot(2,3,4:5)
plot(z-15,amp1(:,jrad,1,1,1)/max(amp1(:,jrad,1,1,1)),'o',z-15,D1CEF(a1(jrad,1,1,1),av/mag,lamem/n,z-15),'--',z-15,D1CEF(a0(jrad,1,1,1),av/mag,lamem/n,z-15));
xlabel('\itz\rm [\mum]')
ylabel('\kappa\rm(\itz\rm)')
