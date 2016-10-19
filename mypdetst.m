function [u,t,x] = mypdetst

m = 0;
x = linspace(-1,1,360);
t = linspace(0,5,20);

sol = pdepe(m,@pdex1pde,@pdex1ic,@pdex1bc,x,t);
% Extract the first solution component as u.
u = sol(:,:,1);

% A surface plot is often a good way to study a solution.
surf(acos(x),t,u)    
title('Numerical solution computed with 100 mesh points.')
xlabel('Angle x')
ylabel('Time t')

% A solution profile can also be illuminating.
figure
plot(acos(x),u)
xlabel('Angle x')
ylabel('u')

% --------------------------------------------------------------
function [c,f,s] = pdex1pde(x,t,u,DuDx)
c = 1;
f = (1-x.^2).*DuDx;
s = -u/(1-x.^2);
% --------------------------------------------------------------
function u0 = pdex1ic(x)
u0 = sqrt(1-x.^2);
% --------------------------------------------------------------
function [pl,ql,pr,qr] = pdex1bc(xl,ul,xr,ur,t)
pl = ul;
ql = 0;
pr = ur;
qr = 0;