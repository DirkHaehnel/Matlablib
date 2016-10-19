function [wp, ep, res, wv] = CylWaveguideMode(nv,rv,m)

% function determines the axial wave vector components for cylindrical waveguide modes
% (c) Joerg Enderlein (2008)

nv = nv(:).'; rv = rv(:).';
dw = 1e-4; % determines the finesse of searching; if waveguide modes are closer than this spacing one has a problem
wv = min(nv)+dw/2:dw:max(nv)-dw/2;

bj = inline('besselj(n,x)','n','x');
bh = inline('besselj(n,x) + i*bessely(n,x)','n','x');
bjx = inline('0.5*(besselj(n-1,x) - besselj(n+1,x))','n','x');
bhx = inline('0.5*(besselj(n-1,x) - besselj(n+1,x) + i*bessely(n-1,x) - i*bessely(n+1,x))','n','x');

len = length(rv);

for jw=1:length(wv)
    w = wv(jw);
    qv = sqrt(nv.^2-w^2);
    qv(imag(qv)<0) = conj(qv(imag(qv)<0));
    if len==1
        r = rv; n2 = nv(2); n1 = nv(1); q2 = qv(2); q1 = qv(1);
        M = [q2.*bhx(m,q2*r), m*w/r.*bh(m,q2*r), -q1.*bjx(m,q1*r), -m*w/r.*bj(m,q1*r); ...
            0, -q2^2*bh(m,q2*r), 0, q1^2*bj(m,q1*r); ...
            i*m*w/r.*bh(m,q2*r), i*n2^2*q2.*bhx(m,q2*r), -i*m*w/r.*bj(m,q1*r), -i*n1^2*q1.*bjx(m,q1*r);
            -i*q2^2*bh(m,q2*r), 0, i*q1^2*bj(m,q1*r), 0 ...
            ];
        %res(jw) = (det([real(M) -imag(M); imag(M) real(M)]));
        res(jw) = real(det(M));
    else
        M = zeros(4*len, 4*len);
        r = rv(end); n2 = nv(end); n1 = nv(end-1); q2 = qv(end); q1 = qv(end-1);
        M(1:4,1:6) = [q2.*bhx(m,q2*r), m*w/r.*bh(m,q2*r), -q1.*bjx(m,q1*r), -m*w/r.*bj(m,q1*r), -q1.*bhx(m,q1*r), -m*w/r.*bh(m,q1*r); ...
            0, -q2^2*bh(m,q2*r), 0, q1^2*bj(m,q1*r), 0, q1^2*bh(m,q1*r); ...
            i*m*w/r.*bh(m,q2*r), i*n2^2*q2.*bhx(m,q2*r), -i*m*w/r.*bj(m,q1*r), -i*n1^2*q1.*bjx(m,q1*r), -i*m*w/r.*bh(m,q1*r), -i*n1^2*q1.*bhx(m,q1*r);
            -i*q2^2*bh(m,q2*r), 0, i*q1^2*bj(m,q1*r), 0, i*q1^2*bh(m,q1*r), 0 ...
            ];
        cnt = 1;
        for k=length(nv)-1:-1:3
            r = rv(k-1); n2 = nv(k); n1 = nv(k-1); q2 = qv(k); q1 = qv(k-1);
            M(1+(4*cnt:4*(cnt+1)),-2+4*cnt+(1:8)) = ...
                [q2.*bjx(m,q2*r), m*w/r.*bj(m,q2*r), q2.*bhx(m,q2*r), m*w/r.*bh(m,q2*r), ...
                -q1.*bjx(m,q1*r), -m*w/r.*bj(m,q1*r), -q1.*bhx(m,q1*r), -m*w/r.*bh(m,q1*r); ...
                0, -q2^2*bj(m,q2*r), 0, -q2^2*bh(m,q2*r), 0, q1^2*bj(m,q1*r), 0, q1^2*bh(m,q1*r); ...
                i*m*w/r.*bj(m,q2*r), i*n2^2*q2.*bjx(m,q2*r), i*m*w/r.*bh(m,q2*r), i*n2^2*q2.*bhx(m,q2*r), ...
                -i*m*w/r.*bj(m,q1*r), -i*n1^2*q1.*bjx(m,q1*r), -i*m*w/r.*bh(m,q1*r), -i*n1^2*q1.*bhx(m,q1*r);
                -i*q2^2*bj(m,q2*r), 0, -i*q2^2*bh(m,q2*r), 0, i*q1^2*bj(m,q1*r), 0, i*q1^2*bh(m,q1*r), 0 ...
                ];
            cnt = cnt +1;
        end
        r = rv(1); n2 = nv(2); n1 = nv(1); q2 = qv(2); q1 = qv(1);
        M(end-3:end,end-5:end) = [q2.*bjx(m,q2*r), m*w/r.*bj(m,q2*r), q2.*bhx(m,q2*r), m*w/r.*bh(m,q2*r), -q1.*bjx(m,q1*r), -m*w/r.*bj(m,q1*r); ...
            0, -q2^2*bj(m,q2*r), 0, -q2^2*bh(m,q2*r), 0, q1^2*bj(m,q1*r); ...
            i*m*w/r.*bj(m,q2*r), i*n2^2*q2.*bjx(m,q2*r), i*m*w/r.*bh(m,q2*r), i*n2^2*q2.*bhx(m,q2*r), -i*m*w/r.*bj(m,q1*r), -i*n1^2*q1.*bjx(m,q1*r);
            -i*q2^2*bj(m,q2*r), 0, -i*q2^2*bh(m,q2*r), 0, i*q1^2*bj(m,q1*r), 0 ...
            ];
        res(jw) = real(det(M));
    end
end

tv = 1:length(wv);
t = tv(res(1:end-1).*res(2:end)<0);
% t = tv(maxl(-res));
wp = []; ep = [];
for j=1:length(t)
    ind = max([1 t(j)-2]):min([length(wv) t(j)+3]);
    %wp(j) = interp1(res(ind).*((tv(ind)<=t(j))-0.5), wv(ind), 0);
    wp(j) = interp1(res(ind), wv(ind), 0, 'cubic');
    q = wp(j);
    qv = sqrt(nv.^2-q^2);
    qv(imag(qv)<0) = conj(qv(imag(qv)<0));
    if len==1
        r = rv; n2 = nv(2); n1 = nv(1); q2 = qv(2); q1 = qv(1);
        M = [q2.*bhx(m,q2*r), m*w/r.*bh(m,q2*r), -q1.*bjx(m,q1*r), -m*w/r.*bj(m,q1*r); ...
            0, -q2^2*bh(m,q2*r), 0, q1^2*bj(m,q1*r); ...
            i*m*w/r.*bh(m,q2*r), i*n2^2*q2.*bhx(m,q2*r), -i*m*w/r.*bj(m,q1*r), -i*n1^2*q1.*bjx(m,q1*r);
            -i*q2^2*bh(m,q2*r), 0, i*q1^2*bj(m,q1*r), 0 ...
            ];
    else
        r = rv(end); n2 = nv(end); n1 = nv(end-1); q2 = qv(end); q1 = qv(end-1);
        M(1:4,1:6) = [q2.*bhx(m,q2*r), m*w/r.*bh(m,q2*r), -q1.*bjx(m,q1*r), -m*w/r.*bj(m,q1*r), -q1.*bhx(m,q1*r), -m*w/r.*bh(m,q1*r); ...
            0, -q2^2*bh(m,q2*r), 0, q1^2*bj(m,q1*r), 0, q1^2*bh(m,q1*r); ...
            i*m*w/r.*bh(m,q2*r), i*n2^2*q2.*bhx(m,q2*r), -i*m*w/r.*bj(m,q1*r), -i*n1^2*q1.*bjx(m,q1*r), -i*m*w/r.*bh(m,q1*r), -i*n1^2*q1.*bhx(m,q1*r);
            -i*q2^2*bh(m,q2*r), 0, i*q1^2*bj(m,q1*r), 0, i*q1^2*bh(m,q1*r), 0 ...
            ];
        cnt = 1;
        for k=length(nv)-1:-1:3
            r = rv(k-1); n2 = nv(k); n1 = nv(k-1); q2 = qv(k); q1 = qv(k-1);
            M(1+(4*cnt:4*(cnt+1)),-2+4*cnt+(1:8)) = ...
                [q2.*bjx(m,q2*r), m*w/r.*bj(m,q2*r), q2.*bhx(m,q2*r), m*w/r.*bh(m,q2*r), ...
                -q1.*bjx(m,q1*r), -m*w/r.*bj(m,q1*r), -q1.*bhx(m,q1*r), -m*w/r.*bh(m,q1*r); ...
                0, -q2^2*bj(m,q2*r), 0, -q2^2*bh(m,q2*r), 0, q1^2*bj(m,q1*r), 0, q1^2*bh(m,q1*r); ...
                i*m*w/r.*bj(m,q2*r), i*n2^2*q2.*bjx(m,q2*r), i*m*w/r.*bh(m,q2*r), i*n2^2*q2.*bhx(m,q2*r), ...
                -i*m*w/r.*bj(m,q1*r), -i*n1^2*q1.*bjx(m,q1*r), -i*m*w/r.*bh(m,q1*r), -i*n1^2*q1.*bhx(m,q1*r);
                -i*q2^2*bj(m,q2*r), 0, -i*q2^2*bh(m,q2*r), 0, i*q1^2*bj(m,q1*r), 0, i*q1^2*bh(m,q1*r), 0 ...
                ];
            cnt = cnt +1;
        end
        r = rv(1); n2 = nv(2); n1 = nv(1); q2 = qv(2); q1 = qv(1);
        M(end-3:end,end-5:end) = [q2.*bjx(m,q2*r), m*w/r.*bj(m,q2*r), q2.*bhx(m,q2*r), m*w/r.*bh(m,q2*r), -q1.*bjx(m,q1*r), -m*w/r.*bj(m,q1*r); ...
            0, -q2^2*bj(m,q2*r), 0, -q2^2*bh(m,q2*r), 0, q1^2*bj(m,q1*r); ...
            i*m*w/r.*bj(m,q2*r), i*n2^2*q2.*bjx(m,q2*r), i*m*w/r.*bh(m,q2*r), i*n2^2*q2.*bhx(m,q2*r), -i*m*w/r.*bj(m,q1*r), -i*n1^2*q1.*bjx(m,q1*r);
            -i*q2^2*bj(m,q2*r), 0, -i*q2^2*bh(m,q2*r), 0, i*q1^2*bj(m,q1*r), 0 ...
            ];
    end
    ep(:,j) = [1; -(M(:,2:end)\M(:,1))];
end


