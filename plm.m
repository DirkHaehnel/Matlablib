function [P,mu,norms]=plm(l,m,mu,check,tol)
% [P,mu,norms]=PLM(l,m,mu,check,tol)
%
% Calculates (associate) Legendre functions, DT (B.48/B.56/B.71).
%
% INPUT:
%
% l      degree (0 <= l <= infinity) [default: random]
% m      order (-l <= m <= l)        [default: all orders 0<=l]
%        l and m can be vectors, but not both at the same time
% mu     argument (-1 <= mu <= 1)    [default: 181 linearly spaced]
% check  1 optional normalization check by Gauss-Legendre quadrature
%        0 no normalization check
% to     Tolerance for optional normalization checking
%
% OUTPUT:
%
% P      The associated Legendre function at the desired argument(s):
%           as a scalar or a row vector with length(mu) columns, OR
%           as a matrix with length(m) rows and length(mu) columns, OR 
%           as a matrix with length(l) rows and length(mu) columns
% mu     The argument(s), which you might or not have specified
% norms  The normalization matrix, which should be the identity matrix
%
% EXAMPLES:
%
% plot(plm([0:5],0,[],1)')
% plot(plm(5,[])')
%
% Last modified by fjsimons-at-alum.mit.edu, 19.05.2006

% Default values
defval('l',round(rand*10))
defval('m',[])
defval('mu',linspace(-1,1,181))
defval('check',0)
defval('tol',1e-10)

% Error handling common to PLM, XLM, YLM
[l,m,mu,check,tol]=pxyerh(l,m,mu,check,tol);

% Warning against using such unnormalized polynomials
% This applies to the PLM's only, really
if 2*l > 21
  warning('Factorials will be inaccurate')
end

switch check
  case 0
   % Calculation for m>=0
   if prod(size(l))==1 & prod(size(m))==1
     % Note that Matlab's unnormalized functions have the (-1)^m
     % Condon-Shortley phase in there so we get rid of it
     P=(-1)^(-m)*rindeks(legendre(l,mu),abs(m)+1);
     if m < 0
       P=(-1)^(-m)*factorial(l+m)./factorial(l-m)*P;
     end
   elseif prod(size(l))==1 
     P=(rindeks(legendre(l,mu),abs(m)+1)'*diag((-1).^(-m)))';
     for index=find(m<0)
       P(index,:)=(-1)^(-m(index))*...
	   factorial(l+m(index))./factorial(l-m(index))*P(index,:);
     end
   elseif prod(size(m))==1
     for index=1:length(l)
       P(index,:)=(-1)^(-m)*rindeks(legendre(l(index),mu),abs(m)+1);
       if m<0
	 P(index,:)=(-1)^(-m)*...
	     factorial(l+m)./factorial(l-m)*P(index,:);
       end
     end
   else
     error('Specify valid option')
   end
   norms=[];
 case 1
  % Check normalization 
  norms=pxynrm(l,m,tol,'P');
  % Still also need to get you the answer at the desired argument
  P=plm(l,m,mu,0);
end
