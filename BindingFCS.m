function [z, zc] = BindingFCS(p, t, y, yc);

% p(1:2) - w^2/4/D in units of delta
% p(3) - delta/sqrt(D);
% p(4) - kp
% p(5) - km

if length(p)==4
    p = p([1 1:4]);
end

t = t(:); 

if p(4)==0
    z = 1./(2*t + sum(p(1:2)));
    zc = z.*exp(-p(3)^2.*z/2);
    z = [z z];
else
    kp = p(4)+p(5);
    km = p(4)-p(5);
    mm = 100;
    q = (0.5:mm)/mm*5/sqrt(p(1));
    for j=1:length(t)
        tmp = sqrt((q.^2+km).^2+4*p(3)*p(4));
        om1 = 0.5*(q.^2+kp-tmp);
        om2 = 0.5*(q.^2+kp+tmp);
        z(j,1) = sum(q.*((kp*(kp+tmp)+km*q.^2).*exp(-om1*t(j))-...
            (kp*(kp-tmp)+km*q.^2).*exp(-om2*t(j)))./(2*kp*tmp).*exp(-q.^2*p(1)));
        z(j,2) = sum(q.*((kp*(kp+tmp)+km*q.^2).*exp(-om1*t(j))-...
            (kp*(kp-tmp)+km*q.^2).*exp(-om2*t(j)))./(2*kp*tmp).*exp(-q.^2*p(2)));
        zc(j) = sum(q.*((kp*(kp+tmp)+km*q.^2).*exp(-om1*t(j))-...
            (kp*(kp-tmp)+km*q.^2).*exp(-om2*t(j)))./(2*kp*tmp).*besselj(0,p(3)*q).*exp(-q.^2*mean(p(1:2))));
    end
end
zc = zc(:)/z(1,1);
z = z/z(1,1);

if nargin>2
    col = ones(length(t),1);
    c1 = [col z(:,1)]\y(:,1); 
    c2 = [col z(:,2)]\y(:,2); 
    zc = [[col z(:,1)]*c1 [col z(:,2)]*c2 [col zc]*sqrt(c1.*c2)];
    semilogx(t,[y yc],'o',t,zc); drawnow;
    z = sum(sum(([y yc]-zc).^2./abs(zc)));
end

return

[z, zc] = BindingFCS([1e-3 1e-2 1e2 1e2], t);
semilogx(t,z,t,1./(t+1e-3)/1e3*z(1))