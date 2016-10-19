% FCSSizeEffect

close all
clear all
a = 0.5; b = 2; % Gaussian parameters of the Gaussian MDF
dx = 1e-2;
[z,rho] = meshgrid(dx/2:dx:3*b,dx/2:dx:3*a);
rad = sqrt(z.^2+rho.^2);

U0 = exp(-2*rho.^2/a^2-2*z.^2/b^2);
[auto0, autotime] = FCS(rho,z,U0,U0);

% generate k-vectors:
maxm = 300;
kv = (0.5:maxm)/maxm*pi/(10*dx);
qv = (0.5:size(U0,2))*pi/max(z(1,:));
% generate bessel and Fourier transform matrices:
Mbessel = besselj(0,kv'*rho(:,1)');
Mfourier = cos(z(1,:)'*qv);
U0wave = Mbessel*(rho.*U0)*Mfourier;

% kk = 5e-5*(kv'.^2*ones(1,length(qv)) + ones(length(kv),1)*qv.^2);
% weight = abs(U0wave).^2.*(kv.'*ones(1,length(qv)));
% for j=1:length(autotime)
%     tst(j) = sum(sum(weight.*exp(-kk*autotime(j))));
% end
% tst = tst/max(tst);

V = 0*rho; V(rad<=0.2) = 1;
Vwave = Mbessel*(rho.*V)*Mfourier;
U1 = Vwave.*U0wave;

U1 = real(Mbessel'*((kv'*ones(size(qv))).*U1)*Mfourier');
auto1 = FCS(rho,z,U1,U1);

V = 0*rho; V(rad>0.19 & rad<=0.21) = 1;
Vwave = Mbessel*(rho.*V)*Mfourier;
U2 = Vwave.*U0wave;
U2 = real(Mbessel'*((kv'*ones(size(qv))).*U2)*Mfourier');
auto2 = FCS(rho,z,U2,U2);

subplot(131)
mpcolor(rho,z,U0,'both'); shading interp; axis image; drawnow
subplot(132)
mpcolor(rho,z,U1,'both'); shading interp; axis image; drawnow
subplot(133)
mpcolor(rho,z,U2,'both'); shading interp; axis image; drawnow


