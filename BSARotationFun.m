function err = BSARotationFun(p, autotime, auto)

beta = sum(auto,5);
t1a = [-flipud(autotime(beta(:,1,1,2)>0)); autotime(beta(:,1,1,1)>0)];
t1b = [-flipud(autotime(beta(:,2,2,2)>0)); autotime(beta(:,2,2,1)>0)];
a1a = [flipud(beta(beta(:,1,1,2)>0,1,1,2)); beta(beta(:,1,1,1)>0,1,1,1)];
a1b = [flipud(beta(beta(:,2,2,2)>0,2,2,2)); beta(beta(:,2,2,1)>0,2,2,1)];
t2 = [-flipud(autotime(beta(:,2,1,2)>0)); autotime(beta(:,1,2,1)>0)];
a2 = [flipud(beta(beta(:,2,1,2)>0,2,1,2)); beta(beta(:,1,2,1)>0,1,2,1)];
t3 = [-flipud(autotime(beta(:,1,2,2)>0)); autotime(beta(:,2,1,1)>0)];
a3 = [flipud(beta(beta(:,1,2,2)>0,1,2,2)); beta(beta(:,2,1,1)>0,2,1,1)];
t = unique([t1a; t1b; t2; t3]);
a1a(abs(t1a)<min(abs(t2))) = [];
a1b(abs(t1b)<min(abs(t2))) = [];
t1a(abs(t1a)<min(abs(t2))) = [];
t1b(abs(t1b)<min(abs(t2))) = [];

if size(autotime,1)<size(auto,1)
else
    t0 = p(1);
    %p = 1./[p(2) 6/20*p(2) p(3)];
    p = 1./[p(2) p(3)];
    z1a = [ones(length(t1a),1) exp(-abs(t1a-t0)*p)];
    z1b = [ones(length(t1b),1) exp(-abs(t1b-t0)*p)];
    z2 = [ones(length(t2),1) exp(-abs(t2-t0)*p)];
    z3 = [ones(length(t3),1) exp(-abs(t3-t0)*p)];
    z = [ones(length(t),1) exp(-abs(t-t0)*p)];
    c1a = z1a\a1a; z1a = z1a*c1a;  zz1a = z*c1a;
    c1b = z1b\a1b; z1b = z1b*c1b; zz1b = z*c1b;
    c2 = z2\a2; z2 = z2*c2; zz2 = z*c2;
    c3 = z3\a3; z3 = z3*c3; zz3 = z*c3;

    [c1a c1b c2 c3]
    
    plot(t1a, a1a, '^r', t1b, a1b, '^r', t2, a2, 'ob', t3, a3, 'vg', t, zz1a, 'y', t, zz1b, 'y', t, zz2, 'c', t, zz3, 'y'); drawnow

    err = sum((a1a-z1a).^2./abs(z1a)) + sum((a1b-z1b).^2./abs(z1b)) + sum((a2-z2).^2./abs(z2)) + sum((a3-z3).^2./abs(z3));
end

