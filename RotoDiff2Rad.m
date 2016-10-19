function rad = RotoDiff2Rad(trot,tt,vv)

if length(tt)<size(trot,2)
    tt = tt*ones(size(trot));
end
tt(tt<100) = tt(tt<100) + 273.15;
T = 273.15 + (0:10:100);   
BoltzmannConstant = 1.380650324e-23; % NIST
if nargin<3 || isempty(vv)
    visc = [1.7920 1.3080 1.0050 0.8007 0.6560 0.5494 0.4688 0.4061 0.3565 0.3165 0.2838];
    vv = interp1(T,visc,tt,'cubic');
end
if size(trot,1)==1
    rad = (6*tt*BoltzmannConstant/8/pi./vv.*trot*1e-6).^(1/3)*1e10; % hydrodyn. radius in Angstrom
else
   da = 1./trot(1,:)/6; % symmetry axis
   db = 1./trot(2,:)/6; % perpendicular axis
   q = 1:0.01:2;
   ft1=sqrt(1-1./q.^2).*q.^(2/3)./log((1+sqrt(1-1./q.^2)).*q); 
   ft2=sqrt(q.^2-1)./q.^(2/3)./atan(sqrt(q.^2-1));
   if da>db
       fra = 4*(1-q.^2)./(3*(2-2*q.^(4/3)./ft2)); % symmetry axis
       frb = 4*(1-q.^4)./q.^2./(3*(2./q.^(2/3).*(2-q.^2)./ft2-2)); % perpendicular axis
   else
       fra = 4*(1-1./q.^2)./(3*(2-2./q.^(4/3)./ft1)); % symmetry axis
       frb = 4*(1-1./q.^4).*q.^2./(3*(2*q.^(2/3).*(2-1./q.^2)./ft1-2)); % perpendicular axis
   end
   fra(1)=1; frb(1)=1;
   q0 = interp1(fra./frb,q,da./db,'cubic');
   tmp = [interp1(q,fra,q0,'cubic'); interp1(q,frb,q0,'cubic')];
   d0 = tmp\[da; db];
   r0 = (6*tt*BoltzmannConstant/8/pi./vv./d0/6*1e-6).^(1/3)*1e10;
   rad = [r0.*q0.^(2/3); r0./q0.^(1/3)];
end
