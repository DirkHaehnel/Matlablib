function [EM,EL,mz,blkm,dblk]=addmout(L)
% [EM,EL,mz,blkm,dblk]=ADDMOUT(L)
%
% Returns spherical harmonic degree and order indexing arrays.
% For arrays where m=-l:l.
%
% INPUT:
%
% L        Maximal degree of the expansion (bandwidth)
%
% OUTPUT:
%
% EM       Vector of orders  m involved in the expansion, [0 -101 -2-1012] 
% EL       Vector of degrees l involved in the expansion, [0 111  222222] 
% mz       Index to the m=0 locations for the zonal coefficients 
% blkm     Reordering sequence for EL/EM to block-sort the orders to 
%          m=0 -1 1 -2 2 ... -L L and then degrees l=abs(m):L within those
% dblk     Reordering sequence to unblock-sort the block-sorted ones...
%
% See also ADDMON, MATRANGES
%
% Last modified by fjsimons-at-alum.mit.edu, 02.05.2006

defval('L',4)

matr=(repmat(0:L,2,1)'*diag([-1 1]))';
EM=matranges(matr(:)')';
twolp=2*(0:L)+1;
EL=gamini(0:L,twolp)';

if nargout>=3
  worist=cumsum(twolp)+1;
  mz=[1 worist(1:end-1)+[0:(L-1)]+1]';
end
if nargout>=4
  % This here is an expression we'd encountered before...
  m=indeks(matr(:)',2:2*(L+1));
  % Find a particular m by the indices:
  matr=[abs(m(:))+1 repmat(L+1,length(m),1)]';
  blkm=[mz(matranges(matr(:)))'+gamini(m,L-abs(m)+1)]';
end
if nargout>=5
  [jk,dblk]=sort(blkm);
end

