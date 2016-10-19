function [coef, M, V] = Cylindrical(qv,mv,nv,rv,c0)

% coef = Cylindrical(qv,mv,nv,rv,c0) calculates the vector-cylindrical
% coefficients coef for scattering the source field c0 at a cylindrical shell
% structure with radii rv and refractive inidces nv
% mv is a vector giving the order of the cylindrical vector functions
% qv is the axial wave vector
%
% order of c0: from out-to-in: m,n/m,n/.../m,n 
% order of coef: from out-to-in: hm,hn/jm,jn,hm,hn/.../jm,jn 

bj = inline('besselj(n,x)','n','x');
bh = inline('besselh(n,x)','n','x');
bjx = inline('0.5*(besselj(n-1,x) - besselj(n+1,x))','n','x');
bhx = inline('0.5*(besselh(n-1,x) - besselh(n+1,x))','n','x');

len = length(rv);
coef = zeros(length(qv),length(mv),4*length(nv)-4);

for k=1:length(qv)
    q = qv(k);
    wv = sqrt(nv.^2-q^2);
    
    if len==1
        r = rv; n2 = nv(2); n1 = nv(1); w2 = wv(2); w1 = wv(1);
        for j=1:length(mv)
            m = mv(j); 
            M = [w2.*bhx(m,w2*r), m*q/r.*bh(m,w2*r), -w1.*bjx(m,w1*r), -m*q/r.*bj(m,w1*r); ...
                    0, -w2.^2.*bh(m,w2*r), 0, w1.^2.*bj(m,w1*r); ...
                    m*q/r.*bh(m,w2*r), n2^2*w2.*bhx(m,w2*r), -m*q/r.*bj(m,w1*r), -n1^2*w1.*bjx(m,w1*r);
                -w2.^2.*bh(m,w2*r), 0, w1.^2.*bj(m,w1*r), 0 ...
                ];
            V = [w2.*bjx(m,w2*r), m*q/r.*bj(m,w2*r), -w1.*bhx(m,w1*r), -m*q/r.*bh(m,w1*r); ...
                    0, -w2.^2.*bj(m,w2*r), 0, w1.^2.*bh(m,w1*r); ...
                    m*q/r.*bj(m,w2*r), n2^2*w2.*bjx(m,w2*r), -m*q/r.*bh(m,w1*r), -n1^2*w1.*bhx(m,w1*r);
                -w2.^2.*bj(m,w2*r), 0, w1.^2.*bh(m,w1*r), 0 ...
                ];
            coef(k,j,:) = (-M\(V*squeeze(c0(k,j,:)))).';
        end
    else
        M = zeros(4*len, 4*len);
        V = M;
        for j=1:length(mv)        
            m = mv(j); 
            r = rv(end); n2 = nv(end); n1 = nv(end-1); w2 = wv(end); w1 = wv(end-1); 
            M(1:4,1:6) = [q.*bhx(m,w2*r), m*w2/r.*bh(m,w2*r), -q.*bjx(m,w1*r), -m*w1/r.*bj(m,w1*r), -q.*bhx(m,w1*r), -m*w1/r.*bh(m,w1*r); ...
                    0, -bh(m,w2*r), 0, bj(m,w1*r), 0, bh(m,w1*r); ...
                    m*w2/r.*bh(m,w2*r), n2^2*q.*bhx(m,w2*r), -m*w1/r.*bj(m,w1*r), -n1^2*q.*bjx(m,w1*r), -m*w1/r.*bh(m,w1*r), -n1^2*q.*bhx(m,w1*r);
                -bh(m,w2*r), 0, bj(m,w1*r), 0, bh(m,w1*r), 0 ...
                ];
            V(1:4,1:6) = [q.*bjx(m,w2*r), m*w2/r.*bj(m,w2*r), -q.*bhx(m,w1*r), -m*w1/r.*bh(m,w1*r), 0, 0; ...
                    0, -bj(m,w2*r), 0, bh(m,w1*r), 0, 0; ...
                    m*w2/r.*bj(m,w2*r), n2^2*q.*bjx(m,w2*r), -m*w1/r.*bh(m,w1*r), -n1^2*q.*bhx(m,w1*r), 0, 0;
                -bj(m,w2*r), 0, bh(m,w1*r), 0, 0, 0 ...
                ];
            cnt = 1;
            for k=length(nv)-1:-1:3
                r = rv(k-1); n2 = nv(k); n1 = nv(k-1); w2 = wv(k); w1 = wv(k-1);        
                M(1+(4*cnt:4*(cnt+1)),-2+4*cnt+(1:8)) = ...
                    [q.*bjx(m,w2*r), m*w2/r.*bj(m,w2*r), q.*bhx(m,w2*r), m*w2/r.*bh(m,w2*r), -q.*bjx(m,w1*r), -m*w1/r.*bj(m,w1*r), -q.*bhx(m,w1*r), -m*w1/r.*bh(m,w1*r); ...
                        0, -bj(m,w2*r), 0, -bh(m,w2*r), 0, bj(m,w1*r), 0, bh(m,w1*r); ...
                        m*w2/r.*bj(m,w2*r), n2^2*q.*bjx(m,w2*r), m*w2/r.*bh(m,w2*r), n2^2*q.*bhx(m,w2*r), -m*w1/r.*bj(m,w1*r), -n1^2*q.*bjx(m,w1*r), -m*w1/r.*bh(m,w1*r), -n1^2*q.*bhx(m,w1*r);
                    -bj(m,w2*r), 0, -bh(m,w2*r), 0, bj(m,w1*r), 0, bh(m,w1*r), 0 ...
                    ];
                V(1+(4*cnt:4*(cnt+1)),-2+4*cnt+(1:8)) = ...
                    [0, 0, q.*bjx(m,w2*r), m*w2/r.*bj(m,w2*r), -q.*bhx(m,w1*r), -m*w1/r.*bh(m,w1*r), 0, 0; ...
                        0, 0, 0, -bj(m,w2*r), 0, bh(m,w1*r), 0, 0; ...
                        0, 0, m*w2/r.*bj(m,w2*r), n2^2*q.*bjx(m,w2*r), -m*w1/r.*bh(m,w1*r), -n1^2*q.*bhx(m,w1*r), 0, 0;
                    0, 0, -bj(m,w2*r), 0, bh(m,w1*r), 0, 0, 0 ...
                    ];
                cnt = cnt +1;
            end
            r = rv(1); n2 = nv(2); n1 = nv(1); w2 = wv(2); w1 = wv(1); 
            M(end-3:end,end-5:end) = [q.*bjx(m,w2*r), m*w2/r.*bj(m,w2*r), q.*bhx(m,w2*r), m*w2/r.*bh(m,w2*r), -q.*bjx(m,w1*r), -m*w1/r.*bj(m,w1*r); ...
                    0, -bj(m,w2*r), 0, -bh(m,w2*r), 0, bj(m,w1*r); ...
                    m*w2/r.*bj(m,w2*r), n2^2*q.*bjx(m,w2*r), m*w2/r.*bh(m,w2*r), n2^2*q.*bhx(m,w2*r), -m*w1/r.*bj(m,w1*r), -n1^2*q.*bjx(m,w1*r);
                -bj(m,w2*r), 0, -bh(m,w2*r), 0, bj(m,w1*r), 0 ...
                ];
            V(end-3:end,end-5:end) = [0, 0, q.*bjx(m,w2*r), m*w2/r.*bj(m,w2*r), -q.*bhx(m,w1*r), -m*w1/r.*bh(m,w1*r); ...
                    0, 0, 0, -bj(m,w2*r), 0, bh(m,w1*r); ...
                    0, 0, m*w2/r.*bj(m,w2*r), n2^2*q.*bjx(m,w2*r), -m*w1/r.*bh(m,w1*r), -n1^2*q.*bhx(m,w1*r);
                0, 0, -bj(m,w2*r), 0, bh(m,w1*r), 0 ...
                ];
            coef(k,j,:) = (-M\(V*squeeze(c0(k,j,:)))).';        
        end
    end
end
