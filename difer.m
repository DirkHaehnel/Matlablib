function difer(V,tolex,sev)
% DIFER(V,tolex,sev)
%
% INPUT:
%
% V         Vector or matrix values
% tolex     Tolerance exponent (default: 10 for 1e-10)
% sev       0 produces WARNING upon failure (default)
%           1 produces ERROR upon failure 
%
% Checks if the sum(abs(V(:))) exceeds 10^(-tolex)
%
% Last modified by fjsimons-at-alum.mit.edu, 18.10.2005

defval('tolex',10)
defval('sev',0)

sabs=sum(abs(V(:)));
if sabs>10^(-tolex)
  mesg=sprintf('sum(abs(%s)) exceeds 0 by %8.3e',inputname(1),sabs);
  switch sev
   case 0
    warning(mesg)
   case 1
    error(mesg)
  end
end



