function val = bint5(n,q,a,b)
% val = bint20(n,q,a,b) integrates the
% function r*besselj(n,q*r) from r=a to 
% r=b using an 5-point Gauss rule.
% If q and a and b are vectors, the result 
% has size [length(max(a,b)) length(q)]

a = a(:)';
b = b(:)';

bp = [-0.90617984593866
  -0.53846931010568
   0.00000000000000
   0.53846931010568
   0.90617984593866];

wf = [0.23692688505619
   0.47862867049937
   0.56888888888889
   0.47862867049937
   0.23692688505619];

val = zeros(length(max(a,b)),length(q));
x = ones(size(bp))*(a+b)/2+bp*(b-a)/2;

for j=1:length(q)
   val(:,j) = (wf'*(x.*besselj(n,q(j)*x).*(ones(size(bp))*(b-a)/2))).'; 
end
