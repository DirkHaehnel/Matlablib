function [err, z] = cornerfun(p, f, y);

f = f(:)';
y = y(:)';
z = 1./(p^2+f.^2);
z = (y/z)*z;

loglog(f,y,'c',f,z); drawnow;

err = sum((y-z).^2./abs(z))
%err = sum((y-z).^2.)