close all

% define interaction functions:
% lambda = 0.5;
% k0 = 2*pi/lambda;
% fun = inline(['cos(' num2str(k0) '*x)./x']);
% fprime = inline(['(cos(' num2str(k0) '*x)./x + ' num2str(k0) '*sin(' num2str(k0) '*x))./x']);
lambda = 0.3;
% fun = inline(['1 + 10*exp(-x/ ' num2str(lambda) '/2) - 11*exp(-x/ ' num2str(lambda) ')']);
% fprime = inline(['10*exp(-x/ ' num2str(lambda) '/2) - 11*exp(-x/ ' num2str(lambda) ')/' num2str(lambda)]);
fun = inline(['exp(-x/' num2str(lambda) ')']);
fprime = inline(['-exp(-x/' num2str(lambda) ')/' num2str(lambda)]);

% generate grid:
if 1
    a = 0.5; b = 2; % Gaussian parameters of the Gaussian MDF
    dx = 1e-2;
    [z,rho] = meshgrid(dx/2:dx:6*b,dx/2:dx:6*a);
    % calculate MDF:
    U0 = exp(-2*rho.^2/a^2 - 2*z.^2/b^2);
else
    rho = exc.rho;
    z = exc.z -exc.focpos(1);
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
kv = (0.5:maxm)/maxm*pi/(10*mean(diff(exc.rho(:,1))));
%qv = (-size(U0,2)+0.5:size(U0,2))*pi/max(z(1,:));
qv = (0.5:size(U0,2)*2)*pi/max(z(1,:))/10;
kk = 5e-2*(kv'.^2*ones(1,length(qv)) + ones(length(kv),1)*qv.^2);

% generate bessel and Fourier transform matrices:
Mbessel = besselj(0,kv'*rho(:,1)');
Mfourier = cos(z(1,:)'*qv);
%Mfourier = exp(-i*z(1,:)'*qv);

U0wave = Mbessel*(rho.*U0)*Mfourier;
% tst = real(Mbessel'*((kv'*ones(size(qv))).*U0wave)*Mfourier');
U0wave2 = abs(U0wave).^2;
cnt = 1;
% V0 = interp1(r0,exp(fun(r0)),rad,'cubic');
% V0(rad<=min(r0))=0;
% nrm = sum(V0(:));
%pcolor(rho,z,U0); shading interp; axis image; drawnow

autotime = [10.^(-2:0.1:4)]; % timescale
for j = 1:length(autotime)
        V = zeros(size(rho));
        V(rad>log(1+autotime(j))) = interp1(r0/2+log(exp(r0)+autotime(j)),exp(-fun(r0)),rad(rad>log(1+autotime(j))),'cubic');
        V(rad<min(0.5*(r0+r))) = 0;
        %V(rad<min(0.5*(r0+r)) | rad>max(0.5*(r0+r))) = 0;
        %V = V/sum(V(:))*nrm;
        Vwave = Mbessel*(rho.*V)*Mfourier;
        %tst = Vwave.*U0wave;
        %tst = real(Mbessel'*((kv'*ones(size(qv))).*tst)*Mfourier');
        %pcolor(rho,z,tst); shading interp; axis image; drawnow
        auto(cnt) = real(sum(sum((kv'*ones(size(qv))).*Vwave.*U0wave2)));
        auto1(cnt) = real(sum(sum((kv'*ones(size(qv))).*Vwave.*U0wave2.*exp(-kk*autotime(j)))));
        auto2(cnt) = real(sum(sum((kv'*ones(size(qv))).*U0wave2.*exp(-kk*autotime(j)))));
        %semilogx(autotime(1:cnt),auto(1:cnt)); drawnow
        surf(V); drawnow
        cnt = cnt + 1;
end

semilogx(autotime,auto/max(auto)); xlabel('time (a.u.)'); ylabel('autocorrelation (norm.)')

%figure; plot(r0,fun(r0)); grid; xlabel('\itr\rm (\mum)'); ylabel('\itV\rm (\itk_B\rm\cdot\itT\rm)')