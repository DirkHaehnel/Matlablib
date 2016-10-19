function err = myerf(x)
% myerf(z) calculates the error function for complex variables

x2=x.*x;
ind = abs(x) < 3.5;
err = x;

if sum(ind(:))>0
    tmp = x2(ind);
    er=1.0;
    r=1.0;
    for  k=1:50;
        r=r.*tmp./(k+0.5);
        er=er+r;
        if prod(double(abs(r) <= abs(er).*eps))==1
            break; 
        end;
    end;
    c0=2.0./sqrt(pi).*x(ind).*exp(-tmp);
    err(ind)=c0.*er;
end

if sum(~ind(:))>0
    tmp = x2(~ind);    
    er=1.0;
    r=1.0;
    for  k=1:12
        r = -r.*(k-0.5)./tmp;
        er=er+r;
    end
    c0 = exp(-tmp)./(abs(x(~ind)).*sqrt(pi));
    tmp = 1.0-c0.*er;
    tmp(real(x(~ind)) < 0.0) = -tmp(real(x(~ind)) < 0.0);
    err(~ind) = tmp;
end;
