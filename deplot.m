function deplot(a,b,l,r,d,t,varargin)

if isempty(r)
    r = l;
end
if nargin<6 || isempty(t)
    t = d;
end
a = a(:)'; b = b(:)'; l = l(:)'; r = r(:)'; d = d(:)'; t = t(:)'; 
for j = 1:length(a)
    xx(1,j) = a(j) - l(j);
    xx(2,j) = a(j) - l(j);
    xx(3,j) = a(j) + r(j);
    xx(4,j) = a(j) + r(j);
    xx(5,j) = a(j) - l(j);
    yy(1,j) = b(j) - d(j);
    yy(2,j) = b(j) + t(j);
    yy(3,j) = b(j) + t(j);
    yy(4,j) = b(j) - d(j);
    yy(5,j) = b(j) - d(j);
end
plot(xx,yy,varargin{:});
