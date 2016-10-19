function Diffraction(dd,nn,pp,ref)

% dd director  

ww = sum(dd.*nn,2);

ts = cross(dd,pp);


tst = abs(ww-1)>1e-5;

dnew(tst,:) = dd(tst,:) + (1/ref-1)*(ww(tst)*[1 1 1]).*nn(tst,:);

dnew = dd./(sqrt(sum(dnew.^2,2))*[1 1 1]);



pp = pp(tst,:) + (1/ref-1)*(sum(pp.*nn,2)*[1 1 1]).*nn(tst,:);

