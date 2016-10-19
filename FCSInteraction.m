function [auto, autotime] = FCSInteraction(p,exc,u)

close all

% define interaction functions:
lambda = p(1);
nukbt = p(2);
fun = inline(['exp(-x/' num2str(lambda) ')']);

% generate grid:
if nargin==1
    a = 0.5; b = 2; % Gaussian parameters of the Gaussian MDF
    dx = 1e-2;
    [z,rho] = meshgrid(dx/2:dx:6*b,dx/2:dx:6*a);
    % calculate MDF:
    U0 = exp(-2*rho.^2/a^2 - 2*z.^2/b^2);
else
    rho = exc.rho;
    z = exc.z - exc.focpos(1);
    U0 = u;
    U0(:,z(1,:)<0) = [];
    rho(:,z(1,:)<0) = [];
    z(:,z(1,:)<0) = [];
end
rad = sqrt(z.^2+rho.^2);
r0 = 0:1e-2:2*max(rad(:))+1e-2; % distances considered for calculating the radial distribution function
r = r0;

% generate k-vectors:
maxm = 300;
kv = (0.5:maxm)/maxm*pi/(10*mean(diff(rho(:,1))));
qv = (0.5:size(U0,2)*2)*pi/max(z(1,:))/10;
kk = 5e-2*(kv'.^2*ones(1,length(qv)) + ones(length(kv),1)*qv.^2);

% generate bessel and Fourier transform matrices:
Mbessel = besselj(0,kv'*rho(:,1)');
Mfourier = cos(z(1,:)'*qv);
%Mfourier = exp(-i*z(1,:)'*qv);

U0wave = Mbessel*(rho.*U0)*Mfourier;
% tst = real(Mbessel'*((kv'*ones(size(qv))).*U0wave)*Mfourier');
U0wave2 = abs(U0wave).^2;
% V0 = interp1(r0,exp(fun(r0)),rad,'cubic');
% V0(rad<=min(r0))=0;
% nrm = sum(V0(:));
%pcolor(rho,z,U0); shading interp; axis image; drawnow

autotime = 10.^(-2:0.1:4); % timescale
auto = zeros(length(autotime),1);
W0 = exp(-nukbt*fun(r0)); W0 = W0/sum(W0);
for j = 1:length(autotime)
        V = zeros(size(rho));
        V(rad>log(1+autotime(j))) = interp1(r0/2+log(exp(r0)+autotime(j)),W0,rad(rad>log(1+autotime(j))),'cubic');
        Vwave = Mbessel*(rho.*V)*Mfourier;
        %tst = Vwave.*U0wave;
        %tst = real(Mbessel'*((kv'*ones(size(qv))).*tst)*Mfourier');
        %pcolor(rho,z,tst); shading interp; axis image; drawnow
        auto(j) = real(sum(sum((kv'*ones(size(qv))).*Vwave.*U0wave2)));
        %auto0(cnt) = real(sum(sum((kv'*ones(size(qv))).*U0wave2.*exp(-kk*autotime(j)))));
        semilogx(autotime(1:j),auto(1:j)); drawnow
        %surf(V); drawnow
end

auto = auto/real(sum(sum((kv'*ones(size(qv))).*U0wave2)));
semilogx(autotime,auto); xlabel('time (a.u.)'); ylabel('autocorrelation (norm.)')

%figure; plot(r0,fun(r0)); grid; xlabel('\itr\rm (\mum)'); ylabel('\itV\rm (\itk_B\rm\cdot\itT\rm)')