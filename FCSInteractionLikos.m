function auto = FCSInteractionLikos(p,autotime,rho,z,u)

u(:,z(1,:)<0) = [];
rho(:,z(1,:)<0) = [];
z(:,z(1,:)<0) = [];

rad = sqrt(z.^2+rho.^2);
dr = rho(2,1)-rho(1,1);
r0 = 10.^(-4:0.01:log10(10*p(1))); % distances considered for calculating the radial distribution function
r = r0;

% generate k-vectors:
maxm = 300;
kv = (0.5:maxm)/maxm*pi/(10*dr);
%qv = (0.5-size(u,2):size(u,2))*pi/max(z(1,:))/10;
qv = (0.5:size(u,2))*pi/max(z(1,:))/10;
kk = 5e-5*(kv'.^2*ones(1,length(qv)) + ones(length(kv),1)*qv.^2);

% generate bessel and Fourier transform matrices:
Mbessel = besselj(0,kv'*rho(:,1)');
Mfourier = cos(z(1,:)'*qv);

uwave = Mbessel*(rho.*u)*Mfourier;
uwave2 = abs(uwave).^2;

fun = inline('exp(-x)./x'); % interaction potential
ei = inline('real(-expint(-x))'); % implicit solution to dynamic equation
eir0 = ei(r0/p(1));
r0(isinf(eir0)) = [];
eir0(isinf(eir0)) = [];

offset = sum(sum(rho.*u.^2));

for j = 1:length(autotime)
    % calculate r(t) as function of r0
    rt = 0*r0;
    ind = eir0+p(2)*autotime(j)<=max(eir0);
    rt(ind) = interp1(eir0(ind),r0(ind),eir0(ind)+p(2)*autotime(j),'cubic');
    rt(~ind) = max(r0);
    
    %plot(r0,rt); pause(0.1);
    
    % calculate R
    rt = (rt+r0)/2;
    
    % calculate r0 corresponding to R
    ind = min(rt(:))<rad & rad<max(rt(:));
    rr = rad;
    rr(ind) = interp1(rt,r0,rad(ind));
    rr(rad<min(rt(:))) = min(rt(:));
    rr(rad>max(rt(:))) = max(rt(:));    
    
    V = exp(-p(3)*fun(rr/p(1)));
    if isempty(min(V(isfinite(V))))
        disp('bla');
    end
    V(~isfinite(V)) = min(V(isfinite(V)));
   
    %surf(V); drawnow

    Vwave = Mbessel*(rho.*V)*Mfourier;
    auto(j) = real(sum(sum((kv'*ones(size(qv))).*Vwave.*uwave2.*exp(-kk*autotime(j)))));
    semilogx(autotime(1:j),auto(1:j)); drawnow
end
auto = auto/max(auto);
semilogx(autotime,auto); xlabel('time (a.u.)'); ylabel('autocorrelation (norm.)')

% plot(r0,fun(r0)); grid; xlabel('\itr\rm (\mum)'); ylabel('\itV\rm (\itk_B\rm\cdot\itT\rm)')