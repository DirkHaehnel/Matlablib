function [err, c, prob] = FIDAFit(p, mn, mh, dn, dh)

% FIDAFit fits the model histogram [mn, mh] against data histogram [dn, dh]

mn = mn(:)'; 
mh = mh(:);
dn = dn(:); dn(dh==0) = [];
dh = dh(:); dh(dh==0) = []; dh = log10(dh);
r = dn(:); rr=r; col = ones(size(r));
rr = cumsum(log([1 1:max(r)])); rr = rr(r+1)'; 
prob = zeros(length(r),length(p));
for j=1:length(p)
    m = round(p(j)*mn); mm=m; mm(m==0) = 1; row = ones(size(m));
    prob(:,j) = real(log10(exp(r*log(mm)-rr*row-col*m)*mh));
end

c = lsqnonneg(prob,dh);

semilogx(dn,dh,dn,prob*c,'o'); drawnow

err = sum((dh - prob*c).^2)

