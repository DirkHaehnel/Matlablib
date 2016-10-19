lambda = 0.5;
n = 1.33;
NA = 1.2;
wd = 3e3;
maxnum = 2e3;
ni = 1.0; % it is assumed that imaging medium is air
ki = 2*pi/lambda*ni; 

% defining ellipsoid of resolution
a = lambda/2/NA;
b = lambda/n/(1-sqrt(1-NA^2/n^2));

chimin = 0;
chimax = asin(NA/n);
dchi = (chimax-chimin)/maxnum;
psi = chimin + (0.5:maxnum)*dchi;
c0 = cos(psi);
s0 = sin(psi);
[v,pc,ps] = DipoleL(psi,2*pi*z/lambda,n1,n,n2,[],2*pi*d/lambda,[]);
v = v.'; ps = ps.'; pc = pc.';

if nargin>12 && ~isempty(atf)
    if length(atf)==1 % accounting for reflection losses when using water immersion; atf = ref index of cover slide
        [~, ~, tmp, tms] = Fresnel(n1(1)*c0,n1(1),atf);
        [~, ~, mp, ms] = Fresnel(sqrt(atf^2-n1(1)^2*s0.^2),atf,n1(1));
    else % + aberration for water immersion; atf(2) = thickness mismatch
        if 1
            [~, ~, tmp, tms] = Fresnel(n1(1)*c0,[n1(1) atf(1) atf(1)], 2*pi*atf(2)/lambda);
            [~, ~, mp, ms] = Fresnel(sqrt(atf(1)^2-n1(1)^2*s0.^2),[atf(1) n1(1) n1(1)], -2*pi*atf(2)/lambda);
        else
            w0 = n1(1)*c0;
            w = sqrt(atf(1)^2-n1(1)^2 + w0.^2);
            mp = exp(1i*2*pi*atf(2)/lambda*(w-w0));
            ms = mp;
            tmp = 1;
            tms = 1;
        end
    end
        v = tmp.*mp.*v;
        pc = tmp.*mp.*pc;
        ps = tms.*ms.*ps;
end

phase = -n1(1).*c0*focpos;
if nargin>13 && ~isempty(ring)
    rad = s0/max(s0);
    eval(['phase = phase + ' ring ';']);
end

%alternative formulation but same result:
%phase = ni*mag*focpos*s0./c0.*si-n1(1)./c0*focpos; 
ez = exp(1i*2*pi/lambda*phase);

fac = dchi*si.*sqrt(ci./c0); % aplanatic objective

barg = ki*si'*rho;
j0 = besselj(0,barg); j1 = besselj(1,barg); j2 = besselj(2,barg);

ezi = fac.*ci.*ez;
ezr = fac.*ez;

fxx0 = (ezi.*pc+ezr.*ps)*j0; % cos(0*phi)-component
fxx2 = -(ezi.*pc-ezr.*ps)*j2; % cos(2*phi)-component
fxz =  -2*1i*(ezi.*v)*j1; % cos(1*phi)-component

byx0 = (ezr.*pc+ezi.*ps)*j0; % cos(0*phi)-component
byx2 = -(ezr.*pc-ezi.*ps)*j2; % cos(2*phi)-component
byz = -2*1i*(ezr.*v)*j1; % cos(1*phi)-component
