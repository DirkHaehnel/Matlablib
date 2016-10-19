function [u,um] = RectDuctVelocity(a,b,y,z,q)
% Velocity profile of rectangular duct given half cross-sectional 
% dimensions (a,b) and coordinates (y,z), and flow (q).

A=0;
for i=1:2:199
A=A+((1/(i.^3))*(-1).^((i-1)/2)*(1-((cosh((i*pi*y)/(2*a)))/(cosh((i*pi*b)/(2*a)))))*cos(i*pi*z/(2*a)));
end
B=0;
for i=1:2:199
B=B+((1/i.^5)*tanh(i*pi*b/(2*a)));
end
area=2*a*2*b;
um=q/area;
c1 = (-3*um)/(((a.^2)*(1-((192/(pi.^5))*(a/b)*B))));
u = ((-16*c1*a.^2)/(pi.^3))*A;