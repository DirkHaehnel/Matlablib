function [ss,N]=gausslegendre(intv,fname,l,method)
% ss=GAUSSLEGENDRE(intv,fname,l,method)
%
% INPUT:
%
% intv       Integration interval
% fname      Inline function name string of integrand
% l          Degree of polynomial accuracy (size 1)
%            or: Matrix with weights and abcissae
% method     Weights calculated by 'cofrec' or 'jacobi'
%
% OUTPUT:
%
% ss         Integral of the integrand
% N          Order of the integration used
%
% For a polynomial function 'fname' and an interval 'intv'
% integrates this function by Gauss-Legendre quadrature.
% You have to give it the polynomial degree of the integrand
% (maximum power, constant is 0), and the integration routine
% will come up with the minimum number of integration points
% for an approximation to the integral which is exact.
% Default will be 10-point integration, but this is pretty
% meaningless. 'FNAME' needs to refer to an existing function
% or to an inline object.
%
% EXAMPLE:
%
% gausslegendre([0 2*pi],'sin')
% gausslegendre([0 1],inline('x.^2'),l);
% gausslegendre([0 pi],inline('x.^2.*cos(x)'),20);
%
% EXAMPLE:
%
% gausslegendre('demo1') % Tests normalization of LEGENDRE
% gausslegendre('demo2') % Tests normalization of LIBBRECHT Schmidt
% gausslegendre('demo3') % Tests normalization of LIBBRECHT Full
% gausslegendre('demo4') % Tests normalization of LIBBRECHT Unnormalized
%
% Last modified by fjsimons-at-alum.mit.edu, April 30th, 2004

if ~isstr(intv)
  defval('l',19)
  defval('method','jacobi')

  %disp('Using Gauss-Legendre integration')

  % Get abcissa's on the interval [-1 1], and weights, and integration
  % order; or receive weights and abcissas as input. If l is a number, need
  % to calculate the weights; if it is a matrix, they were given.
  if length(l)==1
    [w,x,N]=gausslegendrecof(l,method);
  else
    w=l(:,1);
    x=l(:,2);
    N=length(x);
  end

  % Rescale the function values to the interval [a b]
  a=intv(1);
  b=intv(2);
  x=a+(x+1)/2*(b-a);

  % The vectorization is important for degree 0
  fx=feval(fname,x(:)');

  % Calculate integral - note that (b-a)/2 is 1 for [-1 1]
  ss=(w(:)'*fx(:))*(b-a)/2;
else
  if strcmp(intv,'demo1')
  % You will see that the order of integration needs to be exactly
  % the order of the polynomial, or more, and that the accuracy
  % increases in a highly non-linear fashion. At some point, the accuracy
  % decreases again, unless the stable Jacobi algorithm is used!
  deg=round(rand*50);
  g=inline(sprintf('(rindeks(legendre(%i,x),1)).^2',deg));
  lm=100;
  warning off
  for index=1:lm
    gj(index)=gausslegendre([-1 1],g,index-1,'jacobi');
    gc(index)=gausslegendre([-1 1],g,index-1,'cofrec');
  end
  warning on
  plot(0:lm-1,gj-2/(2*deg+1),'+-'); hold on
  plot(0:lm-1,gc-2/(2*deg+1),'ro-'); grid on;
  l=legend('Jacobi','Recursion');
  tl=title(sprintf(...
      'Gauss-Legendre Norm of Unnormalized Legendre function of degree %i',...
      deg));
  xl=xlabel('Degree of polynomial accuracy');
  yl=ylabel('Integration error');
  axis tight
  yli=ylim;
  plot(2*[deg deg],yli,'LineW',2)
  plot(2*[40 40],yli,'LineW',2)
  hold off
  end
  if strcmp(intv,'demo2')
    deg=round(rand*100);
    m=round(rand*deg);
    g=inline(sprintf('(rindeks(libbrecht(%i,x,''sch''),%i)).^2',deg,m+1));
    gj=gausslegendre([-1 1],g,2*deg,'jacobi');
    disp(sprintf(...
      'Norm of Schmidt Legendre function of degree %i and order %i:',...
      deg,m));
    disp(sprintf('%12.10f should be %12.10f; error %8.3e',...
		 gj,(4-2*(m==0))/(2*deg+1),gj-(4-2*(m==0))/(2*deg+1)))
  end
  if strcmp(intv,'demo3')
    deg=round(rand*100);
    m=round(rand*deg);
    g=inline(sprintf('(rindeks(libbrecht(%i,x,''fnc''),%i)).^2',deg,m+1));
    gj=gausslegendre([-1 1],g,2*deg,'jacobi');
    disp(sprintf(...
      'Norm of Full-Norm Legendre function of degree %i and order %i:',...
      deg,m));
    disp(sprintf('%12.10f should be %12.10f; error %8.3e',...
		 gj,1/2/pi,gj-1/2/pi))
  end
end


