function [coef, M, V] = Cylindrical(wv,mv,nv,rv,c0)

% coef = Cylindrical(qv,mv,nv,rv,c0) calculates the vector-cylindrical
% coefficients coef for scattering the source field c0 at a cylindrical shell
% structure with radii rv and refractive inidces nv
% mv is a vector giving the order of the cylindrical vector functions
% qv is the axial wave vector
%
% order of c0: from out-to-in: m,n/m,n/.../m,n 
% order of coef: from out-to-in: hm,hn/jm,jn,hm,hn/.../jm,jn 

warning off MATLAB:nearlySingularMatrix
warning off MATLAB:rankDeficientMatrix

bj = inline('besselj(n,x)','n','x');
%bh = inline('besselj(n,x) + i*bessely(n,x)','n','x');
bh = inline('besselh(n,x)','n','x');
bjx = inline('0.5*(besselj(n-1,x) - besselj(n+1,x))','n','x');
%bhx = inline('0.5*(besselj(n-1,x) - besselj(n+1,x) + i*bessely(n-1,x) - i*bessely(n+1,x))','n','x');
bhx = inline('0.5*(besselh(n-1,x) - besselh(n+1,x))','n','x');

len = length(rv);
coef = zeros(length(wv),length(mv),4*length(nv)-4);

for k=1:length(wv)
    w = wv(k);
    qv = sqrt(nv.^2-w^2);

    if len==1
        r = rv; n2 = nv(2); n1 = nv(1); q2 = qv(2); q1 = qv(1);
        for j=1:length(mv)
            m = mv(j);
            M = [q2.*bhx(m,q2*r), m*w/r.*bh(m,q2*r), -q1.*bjx(m,q1*r), -m*w/r.*bj(m,q1*r); ...
                0, -q2.^2.*bh(m,q2*r), 0, q1.^2.*bj(m,q1*r); ...
                m*w/r.*bh(m,q2*r), n2^2*q2.*bhx(m,q2*r), -m*w/r.*bj(m,q1*r), -n1^2*q1.*bjx(m,q1*r);
                -q2.^2.*bh(m,q2*r), 0, q1.^2.*bj(m,q1*r), 0 ...
                ];
            V = [q2.*bjx(m,q2*r), m*w/r.*bj(m,q2*r), -q1.*bhx(m,q1*r), -m*w/r.*bh(m,q1*r); ...
                0, -q2.^2.*bj(m,q2*r), 0, q1.^2.*bh(m,q1*r); ...
                m*w/r.*bj(m,q2*r), n2^2*q2.*bjx(m,q2*r), -m*w/r.*bh(m,q1*r), -n1^2*q1.*bhx(m,q1*r);
                -q2.^2.*bj(m,q2*r), 0, q1.^2.*bh(m,q1*r), 0 ...
                ];
            coef(k,j,:) = -(M\(V*squeeze(c0(k,j,:)))).';
        end
    else
        M = zeros(4*len, 4*len);
        V = M;
        for j=1:length(mv)        
            m = mv(j); 
            r = rv(end); n2 = nv(end); n1 = nv(end-1); q2 = qv(end); q1 = qv(end-1); 
            M(1:4,1:6) = [q2.*bhx(m,q2*r), m*w/r.*bh(m,q2*r), -q1.*bjx(m,q1*r), -m*w/r.*bj(m,q1*r), -q1.*bhx(m,q1*r), -m*w/r.*bh(m,q1*r); ...
                    0, -q2.^2.*bh(m,q2*r), 0, q1.^2.*bj(m,q1*r), 0, q1.^2.*bh(m,q1*r); ...
                    m*w/r.*bh(m,q2*r), n2^2*q2.*bhx(m,q2*r), -m*w/r.*bj(m,q1*r), -n1^2*q1.*bjx(m,q1*r), -m*w/r.*bh(m,q1*r), -n1^2*q1.*bhx(m,q1*r);
                -q2.^2.*bh(m,q2*r), 0, q1.^2.*bj(m,q1*r), 0, q1.^2.*bh(m,q1*r), 0 ...
                ];
            V(1:4,1:6) = [q2.*bjx(m,q2*r), m*w/r.*bj(m,q2*r), -q1.*bhx(m,q1*r), -m*w/r.*bh(m,q1*r), 0, 0; ...
                    0, -q2.^2.*bj(m,q2*r), 0, q1.^2.*bh(m,q1*r), 0, 0; ...
                    m*w/r.*bj(m,q2*r), n2^2*q2.*bjx(m,q2*r), -m*w/r.*bh(m,q1*r), -n1^2*q1.*bhx(m,q1*r), 0, 0;
                -q2.^2.*bj(m,q2*r), 0, q1.^2.*bh(m,q1*r), 0, 0, 0 ...
                ];
            cnt = 1;
            for kk=length(nv)-1:-1:3
                r = rv(kk-1); n2 = nv(kk); n1 = nv(kk-1); q2 = qv(kk); q1 = qv(kk-1);        
                M(1+4*cnt:4*(cnt+1),-2+4*cnt+(1:8)) = ...
                    [q2.*bjx(m,q2*r), m*w/r.*bj(m,q2*r), q2.*bhx(m,q2*r), m*w/r.*bh(m,q2*r), -q1.*bjx(m,q1*r), -m*w/r.*bj(m,q1*r), -q1.*bhx(m,q1*r), -m*w/r.*bh(m,q1*r); ...
                        0, -q2.^2.*bj(m,q2*r), 0, -q2.^2.*bh(m,q2*r), 0, q1.^2.*bj(m,q1*r), 0, q1.^2.*bh(m,q1*r); ...
                        m*w/r.*bj(m,q2*r), n2^2*q2.*bjx(m,q2*r), m*w/r.*bh(m,q2*r), n2^2*q2.*bhx(m,q2*r), -m*w/r.*bj(m,q1*r), -n1^2*q1.*bjx(m,q1*r), -m*w/r.*bh(m,q1*r), -n1^2*q1.*bhx(m,q1*r);
                    -q2.^2.*bj(m,q2*r), 0, -q2.^2.*bh(m,q2*r), 0, q1.^2.*bj(m,q1*r), 0, q1.^2.*bh(m,q1*r), 0 ...
                    ];
                V(1+4*cnt:4*(cnt+1),-2+4*cnt+(1:8)) = ...
                    [0, 0, q2.*bjx(m,q2*r), m*w/r.*bj(m,q2*r), -q1.*bhx(m,q1*r), -m*w/r.*bh(m,q1*r), 0, 0; ...
                        0, 0, 0, -q2.^2.*bj(m,q2*r), 0, q1.^2.*bh(m,q1*r), 0, 0; ...
                        0, 0, m*w/r.*bj(m,q2*r), n2^2*q2.*bjx(m,q2*r), -m*w/r.*bh(m,q1*r), -n1^2*q1.*bhx(m,q1*r), 0, 0;
                    0, 0, -q2.^2.*bj(m,q2*r), 0, q1.^2.*bh(m,q1*r), 0, 0, 0 ...
                    ];
                cnt = cnt +1;
            end
            r = rv(1); n2 = nv(2); n1 = nv(1); q2 = qv(2); q1 = qv(1); 
            M(end-3:end,end-5:end) = [q2.*bjx(m,q2*r), m*w/r.*bj(m,q2*r), q2.*bhx(m,q2*r), m*w/r.*bh(m,q2*r), -q1.*bjx(m,q1*r), -m*w/r.*bj(m,q1*r); ...
                    0, -q2.^2.*bj(m,q2*r), 0, -q2.^2.*bh(m,q2*r), 0, q1.^2.*bj(m,q1*r); ...
                    m*w/r.*bj(m,q2*r), n2^2*q2.*bjx(m,q2*r), m*w/r.*bh(m,q2*r), n2^2*q2.*bhx(m,q2*r), -m*w/r.*bj(m,q1*r), -n1^2*q1.*bjx(m,q1*r);
                -q2.^2.*bj(m,q2*r), 0, -q2.^2.*bh(m,q2*r), 0, q1.^2.*bj(m,q1*r), 0 ...
                ];
            V(end-3:end,end-5:end) = [0, 0, q2.*bjx(m,q2*r), m*w/r.*bj(m,q2*r), -q1.*bhx(m,q1*r), -m*w/r.*bh(m,q1*r); ...
                    0, 0, 0, -q2.^2.*bj(m,q2*r), 0, q1.^2.*bh(m,q1*r); ...
                    0, 0, m*w/r.*bj(m,q2*r), n2^2*q2.*bjx(m,q2*r), -m*w/r.*bh(m,q1*r), -n1^2*q1.*bhx(m,q1*r);
                0, 0, -q2.^2.*bj(m,q2*r), 0, q1.^2.*bh(m,q1*r), 0 ...
                ];
            tmp = (-M\(V*squeeze(c0(k,j,:)))).';
            if all(isfinite(tmp))
                coef(k,j,:) = tmp;        
            else
                coef(k,j,:) = 0*squeeze(c0(k,j,:)).';        
            end
        end
    end
end

warning on MATLAB:rankDeficientMatrix
warning on MATLAB:nearlySingularMatrix