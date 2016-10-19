function x=indeks(y,i)
% x=INDEKS(y,i)
%
% Extracts indexed positions out of simple matrices
%
% indeks([1 2],2) 
% indeks([1 2],':,2n')
% indeks([1 2],'end')
%
% Works for logical and numeric indices.
%
% Last modified by fjsimons-at-alum.mit.edu, May 3rd, 2004

if ~isstr(i)
  x=y(i);
else
  eval([ 'x=y(' i ');'])
end
