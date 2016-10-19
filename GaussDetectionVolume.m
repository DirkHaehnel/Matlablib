function [v v1 v2] = GaussDetectionVolume(p1, p2, av, lam, sat)

if length(p1)==1
    p1 = [p1 0];
else
    p1 = [p1(1) 0 p1(2)];
end
if length(p2)==2
    p2 = [p2(1) 0 p2(2)];
else
    p2 = [p2 0];
end
w0 = LaserBeamFit(p1,lam(1),0);
if nargin<5 || isempty(sat) || sat==0
    F = @(x) D1CEF(p2,av,lam(2),x);
    v1 = pi*w0^2*quadl(F,0,1e10); 
   
    F = @(x) D1CEF(p2,av,lam(2),x).^2./LaserBeamFit(p1,lam(1),x).^2;
    v2 = pi/2*w0^4*quadl(F,0,1e6);
else
    F = @(x) D1CEF(p2,av,lam(2),x).*log(1+sat*p1(1)^2*D1CEF(p2,av,lam(2),x)./LaserBeamFit(p1,lam(1),x).^2)./...
        (sat*p1(1)^2*D1CEF(p2,av,lam(2),x)./LaserBeamFit(p1,lam(1),x).^2);
    v1 = pi*w0^2*quadl(F,0,1e10); 

    F = @(x) D1CEF(p2,av,lam(2),x).^2./LaserBeamFit(p1,lam(1),x).^2.*...
        ((1+sat*p1(1)^2*D1CEF(p2,av,lam(2),x)./LaserBeamFit(p1,lam(1),x).^2).*...
        log(1+sat*p1(1)^2*D1CEF(p2,av,lam(2),x)./LaserBeamFit(p1,lam(1),x).^2)-...
        sat*p1(1)^2*D1CEF(p2,av,lam(2),x)./LaserBeamFit(p1,lam(1),x).^2)./...
        (sat*p1(1)^2*D1CEF(p2,av,lam(2),x)./LaserBeamFit(p1,lam(1),x).^2).^2./...
        (1+sat*p1(1)^2*D1CEF(p2,av,lam(2),x)./LaserBeamFit(p1,lam(1),x).^2);
    v2 = pi*w0^4*quadl(F,0,1e6);
end
v = v1^2/v2;    
