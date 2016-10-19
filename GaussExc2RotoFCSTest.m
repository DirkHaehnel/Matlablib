function [ux,uy,pp,po,op] = GaussExc2RotoFCSTest

ux = zeros(1,1);
uy = ux; 
pp = zeros(4+1,4+1); 
po = pp; 
op = pp;

[weight,x] = gausslegendrecof(4);
alv = (0:4)/5*2*pi; bev = acos(x);
%clear tst; for j=0:maxL for k=0:j tst(j+1,k+1)=sum(weight'*(SphericalHarmonic(j,k,bev,alv).*conj(SphericalHarmonic(j,k,bev,alv)))); end; end
hwb = waitbar(0,'Calculation in progress');

hxx = 1;
hyy = 0;
hzz = 0;
hxy = 0;
hxz = 0;
hyz = 0;
fxx = 1;
fyy = 0;
fzz = 0;
fxy = 0;
fxz = 0;
fyz = 0;
gyy = 1;
gxx = 0;
gzz = 0;
gxy = 0;
gyz = 0;
gxz = 0;

for jal=1:1;%length(alv)
    al = alv(jal);
    for jbe=1:1;%length(bev)
        be = bev(jbe);
        vx = sin(be)*cos(al);
        vy = sin(be)*sin(al);
        vz = cos(be);
        
        vx = 0;
        vy = 0;
        vz = 1;
        
        poyntx = zeros(5,5,9);
        poynty = poyntx;
        for j=0:0
            for m = 0:j
                for k = -j:j
                    for j1=0:2
                        for m1 = -j1:j1
                            for k1 = -j1:j1
                                for j2=0:2
                                    for m2 = -j2:j2
                                        for k2 = -j2:j2
                                            if abs(j1-j2)<=j && j<=(j1+j2) && m1+m2==m && k1+k2==k
                                                tmp = sqrt((2*j1+1)*(2*j2+1)/(2*j+1))*clebsch(j1,m1,j2,m2,j,m)*clebsch(j1,k1,j2,k2,j,k);
                                                hh = tmp*WignerCoefficients(j2,m2,k2,hxx,hyy,hzz,hxy,hxz,hyz,vx,vy,vz);
                                                poyntx(j+1,m+1,j+1+k) = poyntx(j+1,m+1,j+1+k) + ...
                                                    WignerCoefficients(j1,m1,k1,fxx,fyy,fzz,fxy,fxz,fyz,vx,vy,vz).*hh;
                                                poynty(j+1,m+1,j+1+k) = poynty(j+1,m+1,j+1+k) + ...
                                                    WignerCoefficients(j1,m1,k1,gxx,gyy,gzz,gxy,gxz,gyz,vx,vy,vz).*hh;
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                    %[j m k poyntx(j+1,m+1,j+1+k) poynty(j+1,m+1,j+1+k)]
                end
            end
        end
        ux = ux + weight(jbe)*poyntx(1,1,1);
        uy = uy + weight(jbe)*poynty(1,1,1);
        for L = 0:0
            for M = 0:L
                for K = -L:L
                    pp(L+1,M+1) = pp(L+1,M+1) + (1+sign(M))*weight(jbe)*real(poyntx(L+1,M+1,L+K+1).*conj(poynty(L+1,M+1,L+K+1)));
                    po(L+1,M+1) = po(L+1,M+1) + (1+sign(M))*weight(jbe)*real(poyntx(L+1,M+1,L+K+1).*conj(poyntx(L+1,M+1,L+K+1)));
                    op(L+1,M+1) = op(L+1,M+1) + (1+sign(M))*weight(jbe)*real(poynty(L+1,M+1,L+K+1).*conj(poynty(L+1,M+1,L+K+1)));
                end
            end
        end
        waitbar((jbe+(jal-1)*length(bev))/length(bev)/length(alv),hwb);
    end
end
close(hwb)

% expansion of fxx*vx^2+fyy*vy^2+fzz*vz^2+fxy*vx*vy+fxz*vx*vz+fyz*vy*vz
% in Wigner eigenfunctions of SO(3) group
function f = WignerCoefficients(j,m,k,fxx,fyy,fzz,fxy,fxz,fyz,vx,vy,vz)

switch j
    case 0
        f = (2.*sqrt(2).*(fxx + fyy + fzz).*pi.*(vx.^2 + vy.^2 + vz.^2))/3;
        
    case 1
        f = 0;
        
    case 2
        switch m
            case -2
                switch k
                    case -2
                        f = (fxx + 1i.*fxy - fyy).*pi.*(vx - 1i.*vy).^2/sqrt(10);
                    case -1
                        f = ((fxz + 1i.*fyz).*pi.*(vx - 1i.*vy).^2)/sqrt(10);
                    case 0
                        f = -(((fxx + fyy - 2.*fzz).*pi.*(vx - 1i.*vy).^2)/sqrt(15));
                    case 1
                        f = -(((fxz - 1i.*fyz).*pi.*(vx - 1i.*vy).^2)/sqrt(10));
                    case 2
                        f =  ((fxx - 1i.*fxy - fyy).*pi.*(vx - 1i.*vy).^2)/sqrt(10);
                end
            case -1
                switch k
                    case -2
                        f = sqrt(2/5).*(fxx + 1i.*fxy - fyy).*pi.*(vx - 1i.*vy).*vz;
                    case -1
                        f = sqrt(2/5).*(fxz + 1i.*fyz).*pi.*(vx - 1i.*vy).*vz;
                    case 0
                        f = (-2.*(fxx + fyy - 2.*fzz).*pi.*(vx - 1i.*vy).*vz)/sqrt(15);
                    case 1
                        f = -(sqrt(2/5).*(fxz - 1i.*fyz).*pi.*(vx - 1i.*vy).*vz);
                    case 2
                        f =  sqrt(2/5).*(fxx - 1i.*fxy - fyy).*pi.*(vx - 1i.*vy).*vz;
                end
            case 0
                switch k
                    case -2
                        f = -(((fxx + 1i.*fxy - fyy).*pi.*(vx.^2 + vy.^2 - 2.*vz.^2))/sqrt(15));
                    case -1
                        f = -(((fxz + 1i.*fyz).*pi.*(vx.^2 + vy.^2 - 2.*vz.^2))/sqrt(15));
                    case 0
                        f = (sqrt(2/5).*(fxx + fyy - 2.*fzz).*pi.*(vx.^2 + vy.^2 - 2.*vz.^2))/3;
                    case 1
                        f =  ((fxz - 1i.*fyz).*pi.*(vx.^2 + vy.^2 - 2.*vz.^2))/sqrt(15);
                    case 2
                        f = -(((fxx - 1i.*fxy - fyy).*pi.*(vx.^2 + vy.^2 - 2.*vz.^2))/sqrt(15));
                end
            case 1
                switch k
                    case -2
                        f = -(sqrt(2/5).*(fxx + 1i.*fxy - fyy).*pi.*(vx + 1i.*vy).*vz);
                    case -1
                        f = -(sqrt(2/5).*(fxz + 1i.*fyz).*pi.*(vx + 1i.*vy).*vz);
                    case 0
                        f =  (2.*(fxx + fyy - 2.*fzz).*pi.*(vx + 1i.*vy).*vz)/sqrt(15);
                    case 1
                        f = sqrt(2/5).*(fxz - 1i.*fyz).*pi.*(vx + 1i.*vy).*vz;
                    case 2
                        f = -(sqrt(2/5).*(fxx - 1i.*fxy - fyy).*pi.*(vx + 1i.*vy).*vz);
                end
            case 2
                switch k
                    case -2
                        f = ((fxx + 1i.*fxy - fyy).*pi.*(vx + 1i.*vy).^2)/sqrt(10);
                    case -1
                        f = ((fxz + 1i.*fyz).*pi.*(vx + 1i.*vy).^2)/sqrt(10);
                    case 0
                        f = -(((fxx + fyy - 2.*fzz).*pi.*(vx + 1i.*vy).^2)/sqrt(15));
                    case 1
                        f = -(((fxz - 1i.*fyz).*pi.*(vx + 1i.*vy).^2)/sqrt(10));
                    case 2
                        f = ((fxx - 1i.*fxy - fyy).*pi.*(vx + 1i.*vy).^2)/sqrt(10);
                end
        end
end
