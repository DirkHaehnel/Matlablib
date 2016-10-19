function y = RotoDiffPotentialFun(x,l,m,k,n,j,al,eps)

% RotoDiffPotentialFun helps to calculate the initial 
% multi-vectors for the anisotropic-rotation-in-potential problem
%
% n is the order of the correlator
% j is the number of the correlator element, 1<=j<=n+1
% al is the inclination angle
% eps is the potential energy

if eps>0
    y = binomial(n,j)*eps/(1-exp(-2*eps))*sin(x).*WignerRotation(x,l,m,k).*exp(-eps*(1-cos(x))).*cos(x-al).^(n-j).*sin(x-al).^j;
else
    y = binomial(n,j)/2*sin(x).*WignerRotation(x,l,m,k).*cos(x-al).^(n-j).*sin(x-al).^j;
end
