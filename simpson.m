function egral=simpson(z,f);
% egral=SIMPSON(z,f)
%
% Quadrature on 'f' between z(1) and z(end) using Simpson's rule for
% equally spaced intervals; or trapezoidal rule for unequally spaced
% intervals.
% 'f' may be a matrix here, functions down the columns
%
% Last modified by fjsimons-at-alum.mit.edu, Februrary 9th, 2004

z=z(:);
if size(f,1)==1; f=f(:); end

[m,n]=size(f);

% Odd is the sequence of choice here
% Note that 'z' defines PAIRS
% of layers with equal thickness (as surfc models usually have).
% It checks for that but allows for two variations
% of the equal thickness that then must not exceed
% the next bigger interval/1000.
zd=diff(diff(z)); zd=abs(zd(1:2:end));  zd=zd(~~zd);
if ~isempty(zd)
  uzd=unique(diff(z)); uzd=uzd(3);
  if sum(zd)
    if any(zd>=repmat(uzd/100,length(zd),1))
      disp('Reduced order of integration method to trapezoidal rule')
      egral=trapeze(z,f);
      return
    end
  end
end 

for index=1:n
  egral(index)=[z(3:2:m)-z(1:2:m-2)]'*...
      [f(1:2:m-2,index)+4*f(2:2:m-1,index)+f(3:2:m,index)]/6;
end


    





