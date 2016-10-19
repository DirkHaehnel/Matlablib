function [err, c, z] = GaussFcs(p,NA,av,n,lamem,td,yd,delta)

% p(1) - beam waist
% p(2) - divergence factor

N1 = 1e3;
N2 = 1e2;
tm = 10.^(-3:0.1:5)*p(1)^2;
col = ones(N2,1);
for j=1:length(tm)
    t = tm(j);
    kz = (0:N1-1)/N1*5/sqrt(t);
    z = (0:N1-1)/kz(end)*2*pi;
    q = (0.5:N2)/N2*5/sqrt(t);
    %amp = D3CEF(0,2*pi*z/lamem,2*pi*av/lamem,NA,n); 
    amp = exp(-2*z.^2/5^2/p(1)^2); 
    amp(1) = amp(1)/2;
    ww = p(1)*ones(size(z));%*sqrt(1+(z/p(2)).^2);
    ft = z(2)*fft((col*amp).*exp(-q'.^2*ww.^2/8),[],2);
    ft(1) = ft(1)/sqrt(2);
    ft = kz(2)*abs(ft).^2*exp(-kz'.^2*t);
    ym(j) = diff(q(1:2))*(q.*exp(-q.^2*t))*ft;
end

pd = simplex('AffineFit',1,0,[],[],[],td, yd, tm, ym, 1);
[err, c, z] = AffineFit(pd, td, yd, tm, ym, 1);



