function [fx, fy, fz] = GaussExc2Field(exc,x,y,z)

p = angle(x+i*y);
r = sqrt(x.^2+y.^2);

if nargin>3 && ~isempty(z)
    fx = interp2(exc.rho,exc.z,exc.fxc(:,:,1),r,z,'cubic');
    fy = interp2(exc.rho,exc.z,exc.fyc(:,:,1),r,z,'cubic');
    fz = interp2(exc.rho,exc.z,exc.fzc(:,:,1),r,z,'cubic');
    for j=1:size(exc.fxc,3)-1
        fx = fx + cos(j*p).*interp2(exc.rho,exc.z,exc.fxc(:,:,j+1),r,z,'cubic') + sin(j*p).*interp2(exc.rho,exc.z,exc.fxs(:,:,j),r,z,'cubic');
        fy = fy + cos(j*p).*interp2(exc.rho,exc.z,exc.fyc(:,:,j+1),r,z,'cubic') + sin(j*p).*interp2(exc.rho,exc.z,exc.fys(:,:,j),r,z,'cubic');
        fz = fz + cos(j*p).*interp2(exc.rho,exc.z,exc.fzc(:,:,j+1),r,z,'cubic') + sin(j*p).*interp2(exc.rho,exc.z,exc.fzs(:,:,j),r,z,'cubic');
    end
else
    fx = interp1(exc.rho,exc.fxc(:,1,1),r,'cubic');
    fy = interp1(exc.rho,exc.fyc(:,1,1),r,'cubic');
    fz = interp1(exc.rho,exc.fzc(:,1,1),r,'cubic');
    for j=1:size(exc.fxc,3)-1
        fx = fx + cos(j*p).*interp1(exc.rho,exc.fxc(:,1,j+1),r,'cubic') + sin(j*p).*interp1(exc.rho,exc.fxs(:,1,j),r,'cubic');
        fy = fy + cos(j*p).*interp1(exc.rho,exc.fyc(:,1,j+1),r,'cubic') + sin(j*p).*interp1(exc.rho,exc.fys(:,1,j),r,'cubic');
        fz = fz + cos(j*p).*interp1(exc.rho,exc.fzc(:,1,j+1),r,'cubic') + sin(j*p).*interp1(exc.rho,exc.fzs(:,1,j),r,'cubic');
    end
end
