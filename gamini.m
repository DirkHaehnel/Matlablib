function bigger=gamini(data,folding)
% bigger=GAMINI(data,folding)
%
% INPUT:
%
% data       Some data vector
% folding    The replication factor for every element of the data,
%            if a scalar, this applies to all the elements (default: 3)
%
% EXAMPLE:
%
% [a,b]=degamini(gamini([1 2 3 1 4 5 ],[1 2 3 2 4 2]))
%
% See also GAMINI2
%
% Last modified by fjsimons-at-alum.mit.edu, 15.08.2006

defval('folding',3)

if prod(size(folding))==1
  folding=repmat(folding,prod(size(data)),1);
end

% Never had had a replication factor of 0 before
% But only if they are of the same length!
if length(data)==length(folding)
  data=data(~~folding);
  folding=folding(~~folding);
end

data=data(:)';
folding=folding(:)';

if ~all(size(data)==size(folding)) 
  error([ ' Sizes of input and folding must be the same'])
end

gelp=zeros(1,sum(folding));
gelp([1 cumsum(folding(1:end-1))+1])=1; 
bigger=data(cumsum(gelp));
