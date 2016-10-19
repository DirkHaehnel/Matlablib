function [qp, ep, res, qv] = CylWaveguideMode(nv,rv,m);

% function determines the axial wave vector components for cylindrical waveguide modes
% (c) Joerg Enderlein (2005)

maxval = exp(100);
minval = exp(-100);

nv = nv(:).'; rv = rv(:).';
dq = 1e-4; % determines the finesse of searching; if waveguide modes are closer than this spacing one has a problem
qv = min(nv)+dq/2:dq:max(nv)-dq/2;

bj = inline('besselj(n,x)','n','x');
bh = inline('besselh(n,x)','n','x');
bjx = inline('0.5*(besselj(n-1,x) - besselj(n+1,x))','n','x');
bhx = inline('0.5*(besselh(n-1,x) - besselh(n+1,x))','n','x');

len = length(rv);
coef = [];

for jq=1:length(qv)
    q = qv(jq);
    wv = sqrt(nv.^2-q^2);
    wv(imag(wv)<0) = conj(wv(imag(wv)<0));
    if len==1
        r = rv; n2 = nv(2); n1 = nv(1); w2 = wv(2); w1 = wv(1);
        M = [w2.*bjx(m,w2*r), m*q/r.*bj(m,w2*r), -w1.*bjx(m,w1*r), -m*q/r.*bj(m,w1*r); ...
            0, -w2^2*bh(m,w2*r), 0, w1^2*bj(m,w1*r); ...
            i*m*q/r.*bh(m,w2*r), i*n2^2*w2.*bhx(m,w2*r), -i*m*q/r.*bj(m,w1*r), -i*n1^2*w1.*bjx(m,w1*r);
            -i*w2^2*bh(m,w2*r), 0, i*w1^2*bj(m,w1*r), 0 ...
            ];
        V = [w2.*bhx(m,w2*r), m*q/r.*bh(m,w2*r), -w1.*bhx(m,w1*r), -m*q/r.*bh(m,w1*r); ...
            0, -w2.^2.*bj(m,w2*r), 0, w1.^2.*bh(m,w1*r); ...
            m*q/r.*bj(m,w2*r), n2^2*w2.*bjx(m,w2*r), -m*q/r.*bh(m,w1*r), -n1^2*w1.*bhx(m,w1*r);
            -w2.^2.*bj(m,w2*r), 0, w1.^2.*bh(m,w1*r), 0 ...
            ];

        res(jq) = real(det(M));
        %res(j,jq) = det(M);
    else
        M = zeros(4*len, 4*len);
        r = rv(end); n2 = nv(end); n1 = nv(end-1); w2 = wv(end); w1 = wv(end-1);
        M(1:4,1:6) = [w2.*bhx(m,w2*r), m*q/r.*bh(m,w2*r), -w1.*bjx(m,w1*r), -m*q/r.*bj(m,w1*r), -w1.*bhx(m,w1*r), -m*q/r.*bh(m,w1*r); ...
            0, -w2^2*bh(m,w2*r), 0, w1^2*bj(m,w1*r), 0, w1^2*bh(m,w1*r); ...
            i*m*q/r.*bh(m,w2*r), i*n2^2*w2.*bhx(m,w2*r), -i*m*q/r.*bj(m,w1*r), -i*n1^2*w1.*bjx(m,w1*r), -i*m*q/r.*bh(m,w1*r), -i*n1^2*w1.*bhx(m,w1*r);
            -i*w2^2*bh(m,w2*r), 0, i*w1^2*bj(m,w1*r), 0, i*w1^2*bh(m,w1*r), 0 ...
            ];
        cnt = 1;
        for k=length(nv)-1:-1:3
            r = rv(k-1); n2 = nv(k); n1 = nv(k-1); w2 = wv(k); w1 = wv(k-1);
            M(1+(4*cnt:4*(cnt+1)),-2+4*cnt+(1:8)) = ...
                [w2.*bjx(m,w2*r), m*q/r.*bj(m,w2*r), w2.*bhx(m,w2*r), m*q/r.*bh(m,w2*r), ...
                -w1.*bjx(m,w1*r), -m*q/r.*bj(m,w1*r), -w1.*bhx(m,w1*r), -m*q/r.*bh(m,w1*r); ...
                0, -w2^2*bj(m,w2*r), 0, -w2^2*bh(m,w2*r), 0, w1^2*bj(m,w1*r), 0, w1^2*bh(m,w1*r); ...
                i*m*q/r.*bj(m,w2*r), i*n2^2*w2.*bjx(m,w2*r), i*m*q/r.*bh(m,w2*r), i*n2^2*w2.*bhx(m,w2*r), ...
                -i*m*q/r.*bj(m,w1*r), -i*n1^2*w1.*bjx(m,w1*r), -i*m*q/r.*bh(m,w1*r), -i*n1^2*w1.*bhx(m,w1*r);
                -i*w2^2*bj(m,w2*r), 0, -i*w2^2*bh(m,w2*r), 0, i*w1^2*bj(m,w1*r), 0, i*w1^2*bh(m,w1*r), 0 ...
                ];
            cnt = cnt +1;
        end
        r = rv(1); n2 = nv(2); n1 = nv(1); w2 = wv(2); w1 = wv(1);
        M(end-3:end,end-5:end) = [w2.*bjx(m,w2*r), m*q/r.*bj(m,w2*r), w2.*bhx(m,w2*r), m*q/r.*bh(m,w2*r), -w1.*bjx(m,w1*r), -m*q/r.*bj(m,w1*r); ...
            0, -w2^2*bj(m,w2*r), 0, -w2^2*bh(m,w2*r), 0, w1^2*bj(m,w1*r); ...
            i*m*q/r.*bj(m,w2*r), i*n2^2*w2.*bjx(m,w2*r), i*m*q/r.*bh(m,w2*r), i*n2^2*w2.*bhx(m,w2*r), -i*m*q/r.*bj(m,w1*r), -i*n1^2*w1.*bjx(m,w1*r);
            -i*w2^2*bj(m,w2*r), 0, -i*w2^2*bh(m,w2*r), 0, i*w1^2*bj(m,w1*r), 0 ...
            ];
        res(jq) = real(det(M));
    end
end

t = 1:length(qv);
t = t(res(1:end-1).*res(2:end)<0);
qp = [];
for j=1:length(t)
    ind = max([1 t(j)-2]):min([length(qv) t(j)+3]);
    qp(j) = interp1(res(ind), qv(ind), 0, 'spline');
    q = qp(j);
        wv = sqrt(nv.^2-q^2);
    wv(imag(wv)<0) = conj(wv(imag(wv)<0));
    if len==1
        r = rv; n2 = nv(2); n1 = nv(1); w2 = wv(2); w1 = wv(1);
        M = [w2.*bhx(m,w2*r), m*q/r.*bh(m,w2*r), -w1.*bjx(m,w1*r), -m*q/r.*bj(m,w1*r); ...
            0, -w2^2*bh(m,w2*r), 0, w1^2*bj(m,w1*r); ...
            i*m*q/r.*bh(m,w2*r), i*n2^2*w2.*bhx(m,w2*r), -i*m*q/r.*bj(m,w1*r), -i*n1^2*w1.*bjx(m,w1*r);
            -i*w2^2*bh(m,w2*r), 0, i*w1^2*bj(m,w1*r), 0 ...
            ];
    else
        r = rv(end); n2 = nv(end); n1 = nv(end-1); w2 = wv(end); w1 = wv(end-1);
        M(1:4,1:6) = [w2.*bhx(m,w2*r), m*q/r.*bh(m,w2*r), -w1.*bjx(m,w1*r), -m*q/r.*bj(m,w1*r), -w1.*bhx(m,w1*r), -m*q/r.*bh(m,w1*r); ...
            0, -w2^2*bh(m,w2*r), 0, w1^2*bj(m,w1*r), 0, w1^2*bh(m,w1*r); ...
            i*m*q/r.*bh(m,w2*r), i*n2^2*w2.*bhx(m,w2*r), -i*m*q/r.*bj(m,w1*r), -i*n1^2*w1.*bjx(m,w1*r), -i*m*q/r.*bh(m,w1*r), -i*n1^2*w1.*bhx(m,w1*r);
            -i*w2^2*bh(m,w2*r), 0, i*w1^2*bj(m,w1*r), 0, i*w1^2*bh(m,w1*r), 0 ...
            ];
        cnt = 1;
        for k=length(nv)-1:-1:3
            r = rv(k-1); n2 = nv(k); n1 = nv(k-1); w2 = wv(k); w1 = wv(k-1);
            M(1+(4*cnt:4*(cnt+1)),-2+4*cnt+(1:8)) = ...
                [w2.*bjx(m,w2*r), m*q/r.*bj(m,w2*r), w2.*bhx(m,w2*r), m*q/r.*bh(m,w2*r), ...
                -w1.*bjx(m,w1*r), -m*q/r.*bj(m,w1*r), -w1.*bhx(m,w1*r), -m*q/r.*bh(m,w1*r); ...
                0, -w2^2*bj(m,w2*r), 0, -w2^2*bh(m,w2*r), 0, w1^2*bj(m,w1*r), 0, w1^2*bh(m,w1*r); ...
                i*m*q/r.*bj(m,w2*r), i*n2^2*w2.*bjx(m,w2*r), i*m*q/r.*bh(m,w2*r), i*n2^2*w2.*bhx(m,w2*r), ...
                -i*m*q/r.*bj(m,w1*r), -i*n1^2*w1.*bjx(m,w1*r), -i*m*q/r.*bh(m,w1*r), -i*n1^2*w1.*bhx(m,w1*r);
                -i*w2^2*bj(m,w2*r), 0, -i*w2^2*bh(m,w2*r), 0, i*w1^2*bj(m,w1*r), 0, i*w1^2*bh(m,w1*r), 0 ...
                ];
            cnt = cnt +1;
        end
        r = rv(1); n2 = nv(2); n1 = nv(1); w2 = wv(2); w1 = wv(1);
        M(end-3:end,end-5:end) = [w2.*bjx(m,w2*r), m*q/r.*bj(m,w2*r), w2.*bhx(m,w2*r), m*q/r.*bh(m,w2*r), -w1.*bjx(m,w1*r), -m*q/r.*bj(m,w1*r); ...
            0, -w2^2*bj(m,w2*r), 0, -w2^2*bh(m,w2*r), 0, w1^2*bj(m,w1*r); ...
            i*m*q/r.*bj(m,w2*r), i*n2^2*w2.*bjx(m,w2*r), i*m*q/r.*bh(m,w2*r), i*n2^2*w2.*bhx(m,w2*r), -i*m*q/r.*bj(m,w1*r), -i*n1^2*w1.*bjx(m,w1*r);
            -i*w2^2*bj(m,w2*r), 0, -i*w2^2*bh(m,w2*r), 0, i*w1^2*bj(m,w1*r), 0 ...
            ];
    end
    ep(:,j) = [1; -(M(2:end,2:end)\M(2:end,1))];
end


